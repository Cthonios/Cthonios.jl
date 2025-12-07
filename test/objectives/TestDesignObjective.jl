Random.seed!(123)

# function build_direction_vector(n_design_vars)
function build_direction_vector(params)
    if typeof(params) <: H1Field
        n_design_vars = length(params.data)
    elseif typeof(params) <: NamedTuple
        n_design_vars = mapreduce(length, +, values(params))
    else
        @assert false
    end
    dir_vec = rand(Uniform(-1., 1.), n_design_vars)
    return dir_vec / norm(dir_vec)
end

function generic_dot(a, b)
    if typeof(b) <: H1Field
        return dot(a, b.data)
    elseif typeof(b) <: NamedTuple
        # NOTE this will likely fail for multi-block problems
        # b_vec = reduce(vcat, values(b))
        # return dot(a, b_vec)
        return dot(a, values(b)[1])
    else
        @assert false
    end
end

function perturb_params(params, step_size, direction_vec)
    if typeof(params) <: H1Field
        return H1Field(params + step_size * reshape(direction_vec, size(params)))
    elseif typeof(params) <: NamedTuple
        # return map(
        #     (x, y) -> x .+ step_size * y,
        #     values(params), direction_vec
        # )
        # NOTE this will likely fail for multi-block problems
        # new_params = deepcopy(params)
        # for val in values(new_params)
        #     val .=+ step_size * direction_vec
        # end
        new_params = NamedTuple{keys(params)}((
            values(params)[1] + step_size * direction_vec,
        ))
        return new_params
    else
        @assert false
    end
end

# TODO this is hardcoded for x Sensitivity
function finite_difference_error(
    design_obj, physics_obj, sim,
    step_size, param_sym
)
    objective_cache, U, p = setup_caches(physics_obj, sim; q_degree=1)
    solver = TrustRegionSolver(objective_cache, p; verbose=false)
    des_obj = design_obj(objective_cache, U, p)
    Cthonios.forward_problem!(des_obj, solver, objective_cache, U, p)

    params = getfield(p, param_sym)

    og_obj, og_grad = Cthonios.gradient_and_value(des_obj, params)

    direction_vec = build_direction_vector(params)
    directional_deriv = generic_dot(direction_vec, og_grad)
    params_perturbed = perturb_params(params, step_size, direction_vec)

    objective_cache, U, p = setup_caches(physics_obj, sim; q_degree=1)
    solver = TrustRegionSolver(objective_cache, p; verbose=false)
    des_obj = design_obj(objective_cache, U, p)
    Cthonios.forward_problem!(des_obj, solver, objective_cache, U, p)

    obj_perturbed = Cthonios.value(des_obj, params_perturbed)

    fd_value = (obj_perturbed - og_obj) / step_size
    error = abs(directional_deriv - fd_value)
    return error
end

function strain_energy_objective_with_coords!(f, qois, Us, ps, params)
    U, p = Us[end], ps[end]
    copyto!(p.h1_coords.data, params.data)
    Cthonios._value!(f, qois[1], U, p)
    return nothing
end

function strain_energy_objective_with_props!(f, qois, Us, ps, params)
    U, p = Us[end], ps[end]
    for (n, val) in enumerate(values(params))
        copyto!(values(p.properties)[n], val)
    end
    Cthonios._value!(f, qois[1], U, p)
    return nothing
end

function neohookean_sim_helper()
    mesh_file = dirname(Base.source_dir()) * "/mesh/patch_test_mesh_quad4.g"
    output_file = splitext(mesh_file)[1] * "-output.exo"
    times = TimeStepper(0., 1., 2)
    physics = (;
        var"" = SolidMechanics(
        # Block1 = SolidMechanics(
            PlaneStrain(), NeoHookean()
        )
    )
    props = (;
        var"" = Dict{String, Any}(
        # Block1 = Dict{String, Any}(
            "bulk modulus"  => 1.0,
            "shear modulus" => 0.25
        )
    )

    func_1(x, t) = 0.
    func_2(x, t) = -0.05 * t

    dirichlet_bcs = [
        DirichletBC("displ_x", "sset_1", func_1)
        DirichletBC("displ_y", "sset_1", func_1)
        DirichletBC("displ_x", "sset_3", func_1)
        DirichletBC("displ_y", "sset_3", func_2)
    ]
    objective = Cthonios.QuasiStaticObjective()
    sim = SingleDomainSimulation(
        mesh_file, output_file, 
        times, physics, props;
        dirichlet_bcs=dirichlet_bcs
    )
    return objective, sim
end
 
function neohookean_self_adjoint()
    physics_obj, sim = neohookean_sim_helper()
    qoi = x -> QOIExtractor(
        x, helmholtz_free_energy, +,
        FiniteElementContainers.L2QuadratureField, Float64;
        reduction_2 = +
    )
    design_obj_1 = (x, u, p) -> Cthonios.DesignObjective(
        strain_energy_objective_with_coords!,
        [qoi(x)], u, p
    )
    design_obj_2 = (x, u, p) -> Cthonios.DesignObjective(
        strain_energy_objective_with_props!,
        [qoi(x)], u, p
    )

    for (design_obj, params_sym) in zip(
        [design_obj_1, design_obj_2],
        [:h1_coords, :properties]
    )
        fd_err = finite_difference_error(
            design_obj, physics_obj, sim, 1e-7,
            params_sym
        )
        @test fd_err <= 1e-7

        # TODO setup test for checking v shape
        step_sizes = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
        fd_errs = Float64[]
        for step_size in step_sizes
            fd_err = finite_difference_error(
                design_obj, physics_obj, sim, step_size,
                params_sym
            )
            push!(fd_errs, fd_err)
        end
        display(fd_errs)
    end
end

neohookean_self_adjoint()

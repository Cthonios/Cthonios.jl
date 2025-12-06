Random.seed!(123)

function build_direction_vector(n_design_vars)
    dir_vec = rand(Uniform(-1., 1.), n_design_vars)
    return dir_vec / norm(dir_vec)
end

# TODO this is hardcoded for x Sensitivity
function finite_difference_error(
    forward, 
    objective, sim, qoi,
    step_size
)
    objective_cache, U, p = setup_caches(objective, sim; q_degree=1)
    # sens = forward(objective_cache, U, p, qoi, 2)
    des_obj = forward(objective_cache, U, p, qoi, 2)
    X = p.h1_coords

    # og_obj, og_grad = Cthonios.gradient_x_and_value(sens, U, p)
    # og_obj, og_grad = Cthonios.gradient_and_value(des_obj, U, p)
    og_obj, og_grad = Cthonios.gradient_and_value(des_obj, X)

    direction_vec = build_direction_vector(length(X.data))
    directional_deriv = dot(direction_vec, og_grad.data)

    X_perturbed = H1Field(X + step_size * reshape(direction_vec, size(X)))

    objective_cache, U, p = setup_caches(objective, sim; q_degree=1)
    copyto!(p.h1_coords, X_perturbed)

    # sens = forward(objective_cache, U, p, qoi, 2)
    des_obj = forward(objective_cache, U, p, qoi, 2)
    # obj_perturbed = Cthonios.value(sens, U, p)
    # obj_perturbed = Cthonios.value(des_obj, U, p)
    obj_perturbed = Cthonios.value(des_obj, X_perturbed)

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

function neohookean_forward_problem!(objective_cache, U, p, qoi, n_steps)
    # sens = Cthonios.Sensitivity(qoi(objective_cache), U, p)
    solver = TrustRegionSolver(objective_cache, p)
    Cthonios.initialize!(objective_cache, U, p)

    design_obj = Cthonios.DesignObjective(
        strain_energy_objective_with_coords!,
        [qoi(objective_cache)], U, p
    )

    for _ in 1:n_steps
        Cthonios.step!(solver, objective_cache, U, p)
        push!(design_obj.stored_solutions, deepcopy(U))
        push!(design_obj.stored_parameters, deepcopy(p))
    end

    return design_obj
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
    objective, sim = neohookean_sim_helper()
    qoi = x -> Cthonios.QOIExtractor(
        x, helmholtz_free_energy, +,
        FiniteElementContainers.L2QuadratureField, Float64;
        reduction_2 = +
    )

    
    fd_err = finite_difference_error(
        neohookean_forward_problem!,
        objective, sim, qoi, 1e-7
    )
    @test fd_err <= 1e-7

    # # TODO setup test for checking v shape
    # step_sizes = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
    # fd_errs = Float64[]
    # for step_size in step_sizes
    #     fd_err = finite_difference_error(
    #         neohookean_forward_problem!,
    #         objective, sim, qoi, step_size
    #     )
    #     push!(fd_errs, fd_err)
    # end

    # display(fd_errs)
end

neohookean_self_adjoint()

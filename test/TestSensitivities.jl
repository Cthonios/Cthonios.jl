# using ConstitutiveModels
# using Cthonios
# using Distributions
# using FiniteElementContainers
# using LinearAlgebra
# using Random
# using StaticArrays

Random.seed!(123)

function build_direction_vector(n_design_vars)
    dir_vec = rand(Uniform(-1., 1.), n_design_vars)
    return dir_vec / norm(dir_vec)
end

function finite_difference_error(forward, objective, sim, step_size)
    objective_cache = Cthonios.setup_cache(objective, sim; q_degree=1)
    qoi = Cthonios.QOIExtractor(
        objective_cache, helmholtz_free_energy, +,
        FiniteElementContainers.L2QuadratureField, Float64;
        reduction_2 = +
    )
    sens = forward(objective_cache, qoi, 2)

    U = objective_cache.solution
    p = objective_cache.parameters
    X = p.h1_coords

    og_obj, og_grad = Cthonios.gradient_x_and_value(sens, U, p)

    direction_vec = build_direction_vector(length(X.data))
    directional_deriv = dot(direction_vec, og_grad.data)

    X_perturbed = H1Field(X + step_size * reshape(direction_vec, size(X)))

    objective_cache = Cthonios.setup_cache(objective, sim; q_degree=1)
    copyto!(objective_cache.parameters.h1_coords, X_perturbed)
    qoi = Cthonios.QOIExtractor(
        objective_cache, helmholtz_free_energy, +,
        FiniteElementContainers.L2QuadratureField, Float64;
        reduction_2 = +
    )
    sens = forward(objective_cache, qoi, 2)
    U = objective_cache.solution
    p = objective_cache.parameters
    obj_perturbed = Cthonios.value(sens, U, p)

    fd_value = (obj_perturbed - og_obj) / step_size
    error = abs(directional_deriv - fd_value)
    return error
end

function neohookean_self_adjoint_forward_problem!(objective_cache, qoi, n_steps)
    sens = Cthonios.Sensitivity(qoi)
    solver = TrustRegionSolver(objective_cache)
    Cthonios.initialize!(objective_cache)

    for _ in 1:n_steps
        Cthonios.step!(objective_cache, solver)
        push!(sens.stored_solutions, objective_cache.solution)
        push!(sens.stored_parameters, objective_cache.parameters)
    end

    return sens
end

function neohookean_self_adjoint()
    mesh_file = Base.source_dir() * "/mesh/patch_test_mesh_quad4.g"
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

    sim = SingleDomainSimulation(
        mesh_file, output_file, 
        times, physics, props;
        dirichlet_bcs=dirichlet_bcs
    )
    objective = Cthonios.QuasiStaticObjective()
    fd_err = finite_difference_error(
        neohookean_self_adjoint_forward_problem!,
        objective, sim, 1e-7
    )
    @test fd_err <= 1e-7
end
 
neohookean_self_adjoint()

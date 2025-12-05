using Cthonios
using ConstitutiveModels
using FiniteElementContainers
using LinearAlgebra

# function linear_patch_test_linear_quad_elements()
    mesh_file = Base.source_dir() * "/mesh/path_test_mesh_quad4.g"
    output_file = splitext(mesh_file)[1] * "-output.exo"
    times = TimeStepper(0., 1., 1)
    physics = (;
        var"" = SolidMechanics(
            PlaneStrain(), LinearElastic()
        )
    )
    props = (;
        var"" = Dict{String, Any}(
            "bulk modulus"  => 1.0,
            "shear modulus" => 0.25
        )
    )

    func_1(x, t) = 0.

    dirichlet_bcs = [
        DirichletBC("displ_x", "sset_1", func_1)
        DirichletBC("displ_x", "sset_2", func_1)
        DirichletBC("displ_x", "sset_3", func_1)
        DirichletBC("displ_x", "sset_4", func_1)
        DirichletBC("displ_y", "sset_1", func_1)
        DirichletBC("displ_y", "sset_2", func_1)
        DirichletBC("displ_y", "sset_3", func_1)
        DirichletBC("displ_y", "sset_4", func_1)
    ]

    sim = SingleDomainSimulation(
        mesh_file, output_file, 
        times, physics, props;
        dirichlet_bcs=dirichlet_bcs
    )
    objective = Cthonios.QuasiStaticObjective()
    objective_cache = Cthonios.setup_cache(objective, sim)

    target_disp_grad = [
        0.1 -0.2;
        0.4 -0.1
    ]#' |> collect
    coords = objective_cache.parameters.h1_coords
    # dot(target_disp_grad, coords)
    U = H1Field((coords' * target_disp_grad)' |> collect)
    
    copy!(objective_cache.solution, U)
    solver = Cthonios.TrustRegionSolver(objective_cache)

    Cthonios.solve!(solver, objective_cache.solution.data, objective_cache.parameters)

# end

# linear_patch_test_linear_quad_elements()

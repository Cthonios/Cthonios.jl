function test_quasistatic()
    # file management
    mesh_file = dirname(Base.source_dir()) * "/examples/mechanics/hole_array/mesh/hole_array.exo"
    output_file = splitext(mesh_file)[1] * "-output.exo"

    # Times
    times = TimeStepper(0., 1., 20)

    # Physics
    physics = (;
        Block1 = SolidMechanics(PlaneStrain(), NeoHookean())
    )
    props = (;
        Block1 = Dict{String, Any}(
            "density"         => 1.,
            "Young's modulus" => 1.,
            "Poisson's ratio" => 0.48,
        )
    )

    # Boundary Conditions
    func_1(x, t) = -0.75 * t
    func_2(x, t) = 0.0

    dirichlet_bcs = [
        DirichletBC("displ_x", func_2; sideset_name = "yminus_sideset"),
        DirichletBC("displ_y", func_2; sideset_name = "yminus_sideset"),
        DirichletBC("displ_x", func_2; sideset_name = "yplus_sideset"),
        DirichletBC("displ_y", func_1; sideset_name = "yplus_sideset")
    ]

    # Simulation
    sim = SingleDomainSimulation(
        QuasiStaticObjective,
        mesh_file, output_file, 
        times, physics, props;
        dirichlet_bcs = dirichlet_bcs
    )
    # u = create_unknowns(sim.objective)
    objective = Cthonios.get_objective(sim)
    u, p = Cthonios.get_unknowns_and_parameters(sim)

    timer = TimerOutput()
    preconditioner = CholeskyPreconditioner(objective, u, p, timer)
    predictor = TangentPredictor(objective, u, p, Val{:gmres}(), timer)
    nonlinear_solver = TrustRegionSolver(objective, u, p, timer)
    solver = ImplicitSolver(nonlinear_solver, preconditioner, predictor; timer = timer)
    Cthonios.run!(solver, objective, u, p, sim.mesh, output_file)
end

test_quasistatic()

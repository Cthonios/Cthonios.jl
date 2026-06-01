function test_implicit_dynamics()

    # file management
    mesh_file = Base.source_dir() * "/mesh/sphere.g"
    output_file = splitext(mesh_file)[1] * "-output.exo"

    # Times
    times = TimeStepper(0., 0.01, 11)

    # Physics
    physics = (;
        cube = SolidMechanics(
            ThreeDimensional(), NeoHookean()
        )
    )
    props = (;
        cube = Dict{String, Any}(
            "density"         => 1.e3,
            "Young's modulus" => 1.e4,
            "Poisson's ratio" => 0.33
        )
    )

    # Boundary Conditions
    a = 10
    func_x(x) = -a * x[2] * x[3]
    func_y(x) = a * x[1] * x[3]
    func_z(x) = 0.0

    vel_ics = InitialCondition[
        InitialCondition("displ_x", func_x; nodeset_name = "sphere_surf")
        InitialCondition("displ_y", func_y; nodeset_name = "sphere_surf")
        InitialCondition("displ_z", func_z; nodeset_name = "sphere_surf")
    ]

    # Simulation
    sim = SingleDomainSimulation(
        ImplicitDynamicsObjective,
        mesh_file, output_file, 
        times, physics, props
    )

    # small hack we should remove
    mesh = UnstructuredMesh(mesh_file)
    vel_ics = InitialConditions(mesh, sim.objective.assembler.dof, vel_ics)
    # small hack above we should removes

    Cthonios.initialize!(sim; vel_ics = vel_ics)
    solver = Cthonios.NewtonSolver(sim.objective, sim.p)
    # solver = Cthonios.TrustRegionSolver(objective_cache, p; use_predictor = true)
    Cthonios.run!(sim, solver; output_exodus_every = 10)
end

test_implicit_dynamics()

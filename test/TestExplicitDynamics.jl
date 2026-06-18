function test_explicit_dynamics()
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
        ExplicitDynamicsObjective,
        mesh_file, output_file, 
        times, physics, props
    )
    objective = Cthonios.get_objective(sim)
    u, p = Cthonios.get_unknowns_and_parameters(sim)

    # small hack we should remove
    mesh = UnstructuredMesh(mesh_file)
    vel_ics = InitialConditions(mesh, sim.objective.assembler.dof, vel_ics)
    # small hack above we should removes

    solver = Cthonios.ExplicitSolver(objective, u, p)
    Cthonios.run!(
        solver, objective, u, p, mesh, "test-explicit.exo";
        initialize_settings = (; 
            CFL = 0.1,
            vel_ics = vel_ics
        )
    )
end

test_explicit_dynamics()

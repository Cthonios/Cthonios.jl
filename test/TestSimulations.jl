function test_single_domain_simulation()
    mesh_file = dirname(Base.source_dir()) * "/examples/column_buckle/mesh.g"
    output_file = splitext(mesh_file)[1] * "-output.exo"

    times = TimeStepper(0., 1., 40)

    physics = (;
        block_1 = SolidMechanics(
            ThreeDimensional(), NeoHookean()
        )
    )
    props = (;
        block_1 = Dict{String, Any}(
            "density"       => 1.0,
            "bulk modulus"  => 10.0,
            "shear modulus" => 1.0
        )
    )

    func_1(x, t) = 0.0
    func_2(x, t) = -0.5 * t

    dirichlet_bcs = [
        DirichletBC("displ_x", func_1; sideset_name = "sset_y_negative")
        DirichletBC("displ_y", func_1; sideset_name = "sset_y_negative")
        DirichletBC("displ_z", func_1; sideset_name = "sset_y_negative")
        DirichletBC("displ_x", func_1; sideset_name = "sset_y_positive")
        DirichletBC("displ_z", func_1; sideset_name = "sset_y_positive")
        DirichletBC("displ_y", func_2; sideset_name = "sset_y_positive")
    ]

    sim = SingleDomainSimulation(
        mesh_file, output_file, 
        times, physics, props;
        dirichlet_bcs=dirichlet_bcs
    )

    @test sim.mesh_file == mesh_file
    @test sim.output_file == output_file
    @test sim.times == times
    @test sim.physics == physics
    @test sim.properties == (; block_1 = [1., 10., 1.])
    # TODO add more testing on ics, bcs, and contact pairs
end

test_single_domain_simulation()

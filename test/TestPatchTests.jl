function linear_patch_test_dirichlet(mesh_file, q_degree)
    output_file = splitext(mesh_file)[1] * "-output.exo"
    times = TimeStepper(0., 1., 2)
    physics = (;
        var"" = SolidMechanics(
            PlaneStrain(), LinearElastic()
        )
    )
    props = (;
        var"" = Dict{String, Any}(
            "density"       => 1.0,
            "bulk modulus"  => 1.0,
            "shear modulus" => 0.25
        )
    )

    # NOTE: this bc func is completely incorrect
    # we're hacking things below to conform with
    # the patch test
    func_1(x, t) = 0.
    func_2(x, t) = 0.5 * t

    dirichlet_bcs = [
        DirichletBC("displ_x", func_1; sideset_name = "sset_1")
        DirichletBC("displ_x", func_1; sideset_name = "sset_2")
        DirichletBC("displ_x", func_1; sideset_name = "sset_3")
        DirichletBC("displ_x", func_1; sideset_name = "sset_4")
        DirichletBC("displ_y", func_1; sideset_name = "sset_1")
        DirichletBC("displ_y", func_1; sideset_name = "sset_2")
        DirichletBC("displ_y", func_1; sideset_name = "sset_3")
        DirichletBC("displ_y", func_1; sideset_name = "sset_4")
    ]

    sim = SingleDomainSimulation(
        QuasiStaticObjective,
        mesh_file, output_file, 
        times, physics, props;
        dirichlet_bcs=dirichlet_bcs
    )
    # objective = Cthonios.QuasiStaticObjective()
    # objective_cache, U, p = setup_caches(objective, sim; q_degree=q_degree)

    target_disp_grad = @SMatrix [
        0.1 -0.2;
        0.4 -0.1
    ]
    coords = sim.p.coords
    U = H1Field(target_disp_grad * coords)
    
    # now need to update bcs to reflect what's in U
    # this is hacky but will work for now
    copyto!(sim.p.dirichlet_bcs.bc_cache.vals, U.data[sim.p.dirichlet_bcs.bc_cache.dofs])

    fspace = Cthonios.assembler(sim.objective).dof.var.fspace
    conns = connectivity(fspace, 1)
    ref_fe = FiniteElementContainers.block_reference_element(fspace, 1)

    solver = Cthonios.NewtonSolver(sim.objective, sim.p)

    Uu = create_unknowns(sim.objective)
    FiniteElementContainers.extract_field_unknowns!(Uu, sim.objective.assembler.dof, U)
    Cthonios.solve!(solver, Uu, sim.p)

    grad = Cthonios.gradient(sim.objective, Uu, sim.p)
    @test isapprox(norm(grad), zero(eltype(grad)), atol=1e-12)

    for e in axes(conns, 2)
        for q in 1:ReferenceFiniteElements.num_cell_quadrature_points(ref_fe)
            u_el = FiniteElementContainers._element_level_fields_flat(U, ref_fe, conns, e)
            x_el = FiniteElementContainers._element_level_fields_flat(coords, ref_fe, conns, e)
            interps = FiniteElementContainers._cell_interpolants(ref_fe, q)
            interps = map_interpolants(interps, x_el)
            ∇u_q = interpolate_field_gradients(physics.var"", interps, u_el)
            @test ∇u_q ≈ target_disp_grad
        end
    end
end

linear_patch_test_dirichlet(Base.source_dir() * "/mesh/patch_test_mesh_quad4.g", 1)
linear_patch_test_dirichlet(Base.source_dir() * "/mesh/patch_test_mesh_quad9.g", 2)

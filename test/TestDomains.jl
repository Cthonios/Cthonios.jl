@testset ExtendedTestSet "Domains" begin
  func_1(x, t) = -0.5 * t
  func_2(x, t) = 0.0
  dbcs = [
    DirichletBC("nset_outer_bottom", [1, 2], func_2)
    DirichletBC("nset_outer_top", [1], func_2)
    DirichletBC("nset_outer_top", [2], func_1)
  ]
  nbcs = [

  ]
  sections = Section[
    Section(
      "unnamed_block_1", 2,
      Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
      MaterialProperties(
        "bulk modulus"  => 0.833,
        "shear modulus" => 0.3846
      )
    )
  ]
  domain = Domain("window_pain_tri3.g", sections, dbcs, nbcs)
  time = QuasiStatic(0.0, 1.0, 0.1)
  # test array setup
  coords = domain.coords
  Uu = Cthonios.create_unknowns(domain)
  U = Cthonios.create_fields(domain)
  # TODO add size tests
  @test all(Uu .≈ zero(eltype(Uu)))
  @test all(U .≈ zero(eltype(U)))

  # test assembler setup
  asm = Cthonios.StaticAssembler(domain)
  Cthonios.update_unknown_dofs!(domain, asm)
  ddofs = Int[]
  for bc in domain.dirichlet_bcs
    for dof in bc.dofs
      push!(ddofs, dof)
    end
  end

  for ddof in ddofs
    @test ddof in domain.dirichlet_dofs
  end

  ddofs_test = Cthonios.dirichlet_dofs(domain)
  for ddof in ddofs
    @test ddof in ddofs_test
  end

  Uu = Cthonios.create_unknowns(domain)
  Uu .= rand(eltype(Uu), size(Uu))
  U = Cthonios.create_fields(domain)
  Ubc = zeros(length(ddofs_test))
  Cthonios.step!(time)
  Cthonios.update_dirichlet_vals!(Ubc, domain, coords, time)

  for bc in domain.dirichlet_bcs
    vals = map(x -> bc.func(x, time.current_time[1]), coords[:, bc.nodes])
    for val in vals
      @test val in Ubc
    end
  end

  Cthonios.update_field_bcs!(U, domain, Ubc)
  Cthonios.update_field_unknowns!(U, domain, Uu)
  @test all(U[domain.dirichlet_dofs] .≈ Ubc)
  @test all(U[domain.dof.unknown_dofs] .≈ Uu)
end

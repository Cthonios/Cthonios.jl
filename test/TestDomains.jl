@testset ExtendedTestSet "Domain" begin
  func_1(x, t) = -0.5 * t
  func_2(x, t) = 0.0
  bcs = [
    DirichletBC("nset_outer_bottom", [1, 2], func_2)
    DirichletBC("nset_outer_top", [1], func_2)
    DirichletBC("nset_outer_top", [2], func_1)
  ]
  sections = Section[
    Section(
      Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
      "unnamed_block_1", 2
    )
  ]
  domain = Domain("window_pain_tri3.g", sections, bcs, 2)

  # test array setup
  Uu = Cthonios.create_unknowns(domain)
  U = Cthonios.create_fields(domain)
  # TODO add size tests
  @test all(Uu .≈ zero(eltype(Uu)))
  @test all(U .≈ zero(eltype(U)))

  # test assembler setup
  asm = Cthonios.StaticAssembler(domain)
  Cthonios.update_unknown_dofs!(domain, asm)
end

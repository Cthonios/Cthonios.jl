@testset ExtendedTestSet "Sections" begin
  sections = Section[
    Section(
      Cthonios.Poisson((x, t) -> 1.0),
      "unnamed_block_1", 1
    )
  ]
  @show sections
  @test Cthonios.num_fields(sections[1]) == 1
  @test Cthonios.num_properties(sections[1]) == 0
  @test Cthonios.num_states(sections[1]) == 0

  sections = Section[
    Section(
      Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
      "unnamed_block_1", 2
    )
  ]
  @test Cthonios.num_fields(sections[1].physics) == 2
  @test Cthonios.num_properties(sections[1].physics) == 2
  @test Cthonios.num_states(sections[1].physics) == 0
  @test sections[1].physics.material_model == NeoHookean()
  @test sections[1].physics.formulation == PlaneStrain()

  # sections = map()
  bcs = [
    DirichletBC("nset_outer_bottom", [1, 2], (x, t) -> 0.0),
    DirichletBC("nset_outer_top", [1], (x, t) -> 0.0),
    DirichletBC("nset_outer_top", [2], (x, t) -> -0.5 * t)
  ]
  mesh = FileMesh(ExodusDatabase, "window_pain_tri3.g")
  bcs = map(bc -> Cthonios.DirichletBCInternal(mesh, bc, 2), bcs)
  coords = coordinates(mesh)
  coords = NodalField{size(coords), Vector}(coords)
  dof = DofManager{2, size(coords, 2), Vector{Float64}}()

  section = Cthonios.SectionInternal(mesh, dof, sections[1])
  @show section
  @test Cthonios.num_fields(section) == 2
  @test Cthonios.num_properties(section) == 2
  @test Cthonios.num_states(section) == 0
  @test section.physics.material_model == NeoHookean()
  @test section.physics.formulation == PlaneStrain()
  @test size(section) == (2, 3, 2, 0)
end
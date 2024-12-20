@testset ExtendedTestSet "Dirichlet Boundary Conditions" begin
  bc_input = DirichletBC("nset_inner_top", [1], (x, t) -> 0.0)
  @show bc_input
  @test bc_input.nset_name == "nset_inner_top"
  @test bc_input.dofs == [1]
  @test bc_input.func(SVector{2, Float64}(1., 2.), 1.) == 0.0

  bc_input = DirichletBC("nset_inner_top", [1, 2], (x, t) -> 1.0)
  @test bc_input.nset_name == "nset_inner_top"
  @test bc_input.dofs == [1, 2]
  @test bc_input.func(SVector{2, Float64}(1., 2.), 1.) == 1.0



  mesh = FileMesh(ExodusDatabase, "window_pain_tri3.g")
  bc = Cthonios.DirichletBCInternal(mesh, bc_input, 2)
  @show bc

  # bad input
  bc_input = DirichletBC("nset_bad", [1, 2], (x, t) -> 1.0)
  @test_throws Exodus.SetNameException bc = Cthonios.DirichletBCInternal(mesh, bc_input, 2)

  bc_input = DirichletBC("nset_bad", [1, 2, 3], (x, t) -> 1.0)
  @test_throws ErrorException bc = Cthonios.DirichletBCInternal(mesh, bc_input, 2)
end

@testset ExtendedTestSet "Neumann Boundary Conditions" begin
  bc_input = NeumannBC("sset_inner_top", (x, t) -> (0., 0.))
  @show bc_input
  @test bc_input.sset_name == "sset_inner_top"
  @test bc_input.func(SVector{2, Float64}(1., 2.), 1.) == (0., 0.)

  mesh = FileMesh(ExodusDatabase, "window_pain_tri3.g")
  bc = Cthonios.NeumannBCInternal(mesh, bc_input)
  @show bc

  # bad input
  bc_input = NeumannBC("sset_bad", (x, t) -> (1.0, 0.))
  @test_throws Exodus.SetNameException bc = Cthonios.NeumannBCInternal(mesh, bc_input)
end
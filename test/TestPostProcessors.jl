@testset ExtendedTestSet "PostProcessors" begin
  mesh = FileMesh(ExodusDatabase, "window_pain_tri3.g")
  coords = coordinates(mesh)
  coords = NodalField{size(coords), Vector}(coords)
  dof = DofManager{2, size(coords, 2), Vector{Float64}}()
  U1 = FiniteElementContainers.create_fields(dof)
  U2 = FiniteElementContainers.create_fields(dof)
  U1 .= rand(eltype(U1), size(U1))
  U2 .= rand(eltype(U2), size(U2))


  pp = ExodusPostProcessor("window_pain_tri3.g", "output.e", ["displ_x", "displ_y"])
  Cthonios.write_time(pp, 1, 0.0)
  Cthonios.write_time(pp, 2, 1.0)
  @test Exodus.read_times(pp.exo)[1] ≈ 0.0
  @test Exodus.read_times(pp.exo)[2] ≈ 1.0
  Cthonios.write_fields(pp, U1, 1)
  Cthonios.write_fields(pp, U2, 2)
  @test all(Exodus.read_values(pp.exo, NodalVariable, 1, "displ_x") .≈ U1[1, :])
  @test all(Exodus.read_values(pp.exo, NodalVariable, 1, "displ_y") .≈ U1[2, :])
  @test all(Exodus.read_values(pp.exo, NodalVariable, 2, "displ_x") .≈ U2[1, :])
  @test all(Exodus.read_values(pp.exo, NodalVariable, 2, "displ_y") .≈ U2[2, :])
  Cthonios.close(pp)
end

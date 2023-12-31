using Aqua
using Cthonios
using JET
using Test
using TestSetExtensions

@testset "Cthonios.jl" begin
    # Write your tests here.
end

@testset ExtendedTestSet "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false)
end

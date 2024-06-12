using Aqua
using Cthonios
using JET
using Test
using TestSetExtensions

@testset ExtendedTestSet "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false)
end

@testset ExtendedTestSet "JET.jl" begin
  JET.test_package(Cthonios; target_defined_modules=true)
end

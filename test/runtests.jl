using Aqua
using Cthonios
using Exodus
using JET
using Test
using TestSetExtensions

# for regression_test in readdir("regression_tests/")
#   @testset ExtendedTestSet "$regression_test" begin
#     dir = joinpath("regression_tests/", regression_test)
#     Cthonios.cthonios_main(joinpath(dir, "test.yaml"), false)
#     @test Exodus.exodiff(joinpath(dir, "test.e"), joinpath(dir, "test.gold.e"))
#     rm(joinpath(dir, "test.e"); force=true)
#     rm(joinpath(dir, "test.log"); force=true)
#     rm("exodiff.log"; force=true)
#   end
# end

@testset ExtendedTestSet "Newton - Tri3" begin
  Cthonios.cthonios_main("regression_tests/newton_tri3_uniaxial_tension/test.yaml", false, "NoKABackend")
end

@testset ExtendedTestSet "Trust region - Tri3" begin
  Cthonios.cthonios_main("regression_tests/trs_tri3_uniaxial_tension/test.yaml", false, "NoKABackend")
end

@testset ExtendedTestSet "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false)
end

@testset ExtendedTestSet "JET.jl" begin
  JET.test_package(Cthonios; target_defined_modules=true)
end

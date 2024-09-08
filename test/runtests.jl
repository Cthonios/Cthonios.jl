using Aqua
using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using JET
using StaticArrays
using Test
using TestSetExtensions

include("TestBoundaryConditions.jl")
include("TestDomains.jl")
include("TestPostProcessors.jl")
include("TestSections.jl")
include("TestTimeSteppers.jl")

@testset ExtendedTestSet "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false)
end

@testset ExtendedTestSet "JET.jl" begin
  JET.test_package(Cthonios; target_defined_modules=true)
end

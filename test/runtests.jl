using Aqua
using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using JET
using StaticArrays
using Test
using TestSetExtensions
using TimerOutputs

include("TestBoundaryConditions.jl")
include("TestDomains.jl")
include("TestIntegrators.jl")
include("TestObjectives.jl")
include("TestPostProcessors.jl")
include("TestProblems.jl")
include("TestSections.jl")

@testset ExtendedTestSet "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false, persistent_tasks=false)
end

@testset ExtendedTestSet "JET.jl" begin
  JET.test_package(Cthonios; target_defined_modules=true)
end

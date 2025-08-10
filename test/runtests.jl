using Aqua
using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using StaticArrays
using Test
using TestSetExtensions
using TimerOutputs

# include("TestBoundaryConditions.jl")
# include("TestDomains.jl")
# include("TestIntegrators.jl")
# include("TestObjectives.jl")
# include("TestPostProcessors.jl")
# include("TestProblems.jl")
# include("TestSections.jl")

@testset ExtendedTestSet "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false, persistent_tasks=false)
end

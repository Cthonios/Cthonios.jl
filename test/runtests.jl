using Aqua
using ConstitutiveModels
using Cthonios
using FiniteElementContainers
using ForwardDiff
using LinearAlgebra
using StaticArrays
using Test

@testset "CLI" begin
  include("TestCLI.jl")
end

@testset "Contact" begin
  include("TestContact.jl")
end

@testset "Objectives" begin
  include("TestObjectives.jl")
end

@testset "Physics" begin
  include("TestPhysics.jl")
end

@testset "PostProcessors" begin
  include("TestPostProcessors.jl")
end

@testset "Simulations" begin
  include("TestSimulations.jl")
end

@testset "Solvers" begin
  include("TestSolvers.jl")
end

@testset "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false, persistent_tasks=false)
end

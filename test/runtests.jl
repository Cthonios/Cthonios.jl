using Aqua
using ConstitutiveModels
using Cthonios
using FiniteElementContainers
using StaticArrays
using Test

@testset "Objectives" begin
  include("TestObjectives.jl")
end

@testset "Physics" begin
  include("TestPhysics.jl")
end

@testset "Simulations" begin
  include("TestSimulations.jl")
end

@testset "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false, persistent_tasks=false)
end

using Aqua
using ConstitutiveModels
using Cthonios
using FiniteElementContainers
using StaticArrays
using Test

@testset "Physics.jl" begin
  include("TestPhysics.jl")
end

@testset "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false, persistent_tasks=false)
end

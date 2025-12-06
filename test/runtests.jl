using Aqua
using ConstitutiveModels
using Cthonios
using Distributions
using FiniteElementContainers
using ForwardDiff
using LinearAlgebra
using Random
using StaticArrays
using Test

function sim_helper()
  mesh_file = dirname(dirname(@__FILE__)) * "/examples/column_buckle/mesh.g"
  output_file = splitext(mesh_file)[1] * "-output.exo"

  times = TimeStepper(0., 1., 40)

  physics = (;
      block_1 = SolidMechanics(
          ThreeDimensional(), NeoHookean()
      )
  )
  props = (;
      block_1 = Dict{String, Any}(
      "bulk modulus"  => 10.0,
      "shear modulus" => 1.0
      )
  )

  func_1(x, t) = 0.0
  func_2(x, t) = -0.5 * t

  dirichlet_bcs = [
      DirichletBC("displ_x", "sset_y_negative", func_1)
      DirichletBC("displ_y", "sset_y_negative", func_1)
      DirichletBC("displ_z", "sset_y_negative", func_1)
      DirichletBC("displ_x", "sset_y_positive", func_1)
      DirichletBC("displ_z", "sset_y_positive", func_1)
      DirichletBC("displ_y", "sset_y_positive", func_2)
  ]

  sim = SingleDomainSimulation(
      mesh_file, output_file, 
      times, physics, props;
      dirichlet_bcs=dirichlet_bcs
  )
  return sim
end

@testset "CLI" begin
  include("TestCLI.jl")
end

@testset "Contact" begin
  include("TestContact.jl")
end

@testset "Objectives" begin
  include("objectives/TestDesignObjective.jl")
  include("objectives/TestObjectives.jl")
end

@testset "Patch Tests" begin
  include("TestPatchTests.jl")
end

@testset "Physics" begin
  include("TestPhysics.jl")
end

@testset "PostProcessors" begin
  include("TestPostProcessors.jl")
end

# @testset "Sensitivities" begin
#   include("TestSensitivities.jl")
# end

@testset "Simulations" begin
  include("TestSimulations.jl")
end

@testset "Solvers" begin
  include("TestSolvers.jl")
end

@testset "Aqua.jl" begin
  Aqua.test_all(Cthonios; ambiguities=false, persistent_tasks=false)
end

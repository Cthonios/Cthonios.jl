module Cthonios

export DirectSolver
export DofManager, ExodusDatabase, FileMesh
export MatrixFreeSolver, NewtonSolver
export Objective
export PlaneStrain
export PostProcessor
export TotalLagrangeSection
# export create_unknowns

export DirichletBC, NeumnanBC
export Domain
export TotalLagrangeSection
export ConstantTimeStepper

import Exodus: ExodusDatabase
import FiniteElementContainers: DofManager, FileMesh, PlaneStrain

using ComponentArrays
using ConstitutiveModels
using DocStringExtensions
using Exodus
using FiniteElementContainers
using LinearAlgebra
using LinearMaps
using LinearSolve
using Printf
using RuntimeGeneratedFunctions
using SparseArrays
using StaticArrays
using TimerOutputs
using UnPack
using YAML

RuntimeGeneratedFunctions.init(@__MODULE__)

function add_tabs(str, num_tabs)
  if num_tabs == 0
    return str
  end
  lines = split(repr(str), "\n")
  for n in axes(lines, 1)
    lines[n] = repeat("  ", num_tabs) * lines[n] * "\n"
  end
  # @show lines
  return string(lines...)
end

function input_with_default(inputs::Dict{Symbol, Any}, key::String, default)
  key = Symbol(key)
  if key in keys(inputs)
    val = inputs[key]
  else
    val = default
  end
  val
end

include("BoundaryConditions.jl")
include("Parsers.jl")
include("PostProcessors.jl")
include("Sections.jl")
include("TimeSteppers.jl")

include("Domains.jl")
include("Mechanics.jl")
include("Solvers.jl")

include("Objective.jl")
include("Problems.jl")

end # module

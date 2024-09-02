module Cthonios

using ComponentArrays
using ConstitutiveModels
using DocStringExtensions
using Exodus
using FiniteElementContainers
using LinearAlgebra
using LinearSolve
using Parameters
using Printf
using RuntimeGeneratedFunctions
using SparseArrays
using StaticArrays
using TimerOutputs

RuntimeGeneratedFunctions.init(@__MODULE__)

include("BoundaryConditions.jl")
include("physics/Physics.jl")
include("PostProcessors.jl")
include("Sections.jl")
include("TimeSteppers.jl")

include("Domains.jl")
include("Objectives.jl")

include("Iterators.jl")
include("solvers/Solvers.jl")

include("Problems.jl")

end # module

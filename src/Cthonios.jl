module Cthonios

using ArgParse
using Atomix
using ComponentArrays
using ConstitutiveModels
using DifferentiationInterface
using DocStringExtensions
using Exodus
using FiniteElementContainers
using IterativeSolvers
using LimitedLDLFactorizations
using LDLFactorizations
using LinearAlgebra
using LinearSolve
using Parameters
using Printf
using ReferenceFiniteElements
using RuntimeGeneratedFunctions
using SciMLOperators
using SparseArrays
using StaticArrays
using TimerOutputs
using YAML

RuntimeGeneratedFunctions.init(@__MODULE__)

include("bcs/BoundaryConditions.jl")
include("physics/Physics.jl")
include("PostProcessors.jl")
include("Sections.jl")
include("TimeSteppers.jl")

include("Domains.jl")

include("contact/Contact.jl")

include("Objectives.jl")

include("iterators/Iterators.jl")
include("solvers/Solvers.jl")

include("problems/Problems.jl")

include("Main.jl")

end # module

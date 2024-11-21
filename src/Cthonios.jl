module Cthonios

using Reexport

using ArgParse
using Atomix
using ComponentArrays
@reexport using ConstitutiveModels
using DocStringExtensions
@reexport using Exodus
@reexport using FiniteElementContainers
using ForwardDiff
using IterativeSolvers
using LimitedLDLFactorizations
using LDLFactorizations
@reexport using LinearAlgebra
using LinearSolve
using Parameters
using Printf
using ReferenceFiniteElements
using RuntimeGeneratedFunctions
using SciMLOperators
using SparseArrays
@reexport using StaticArrays
@reexport using TimerOutputs
using YAML

# import to avoid name conflicts with gradient, hvp
import Enzyme: Const, Duplicated, Forward, Reverse, autodiff, make_zero, make_zero!

RuntimeGeneratedFunctions.init(@__MODULE__)

include("bcs/BoundaryConditions.jl")
include("physics/Physics.jl")
include("PostProcessors.jl")
include("Sections.jl")
# include("TimeSteppers.jl")

include("Domains.jl")

include("contact/Contact.jl")

include("integrators/Integrators.jl")
include("Objectives.jl")

include("iterators/Iterators.jl")
include("solvers/Solvers.jl")

# include("problems/Problems.jl")
include("Problem.jl")

include("Main.jl")

end # module

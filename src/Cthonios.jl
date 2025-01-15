module Cthonios

# import this first since we're using it below
using Reexport

using ArgParse
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
using Parameters
using Printf
using ReferenceFiniteElements
using RuntimeGeneratedFunctions
using SparseArrays
@reexport using StaticArrays
@reexport using TimerOutputs
using YAML

# import to avoid name conflicts with gradient, hvp
import Enzyme
import Enzyme: Const, Duplicated, Forward, Reverse, autodiff, make_zero, make_zero!

RuntimeGeneratedFunctions.init(@__MODULE__)

# small components that make up a more complex problem
include("bcs/BoundaryConditions.jl")
include("physics/Physics.jl")
include("PostProcessors.jl")
include("Sections.jl")

# container for small components
include("Domains.jl")

# extra optional stuff related to constraints
include("contact/Contact.jl")

# time integration facilities
include("integrators/Integrators.jl")

# integrals over objective function kernels
include("integrals/Integrals.jl")

# different objectives
include("objectives/Objectives.jl")

# integrals over objective function kernels
# include("integrals/Integrals.jl")

# solvers
include("solvers/Solvers.jl")

# different problem types
include("Problem.jl")
# include("SchwarzProblem.jl")

# CLI
include("Main.jl")

end # module

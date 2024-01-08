module Cthonios

# edpxort
export CthoniosCommon

# export create_field

# dev stuff

# export TotalLagrangeSection


# dependencies
using ArgParse
using ConstitutiveModels
using DocStringExtensions
using Exodus
using FiniteElementContainers
using FunctionWrappers
using IterativeSolvers
using LinearAlgebra
using LinearSolve
using Logging
using LoggingExtras
using Parameters
using Preconditioners
using Printf
using ReferenceFiniteElements
using StaticArrays
using StructArrays
using Tensors
using TimerOutputs
using YAML

import FunctionWrappers: FunctionWrapper


# for docs
@template (FUNCTIONS, METHODS, MACROS) = 
"""
$(TYPEDSIGNATURES)
$(DOCSTRING)
$(METHODLIST)
"""

@template (TYPES) = 
"""
$(TYPEDFIELDS)
$(DOCSTRING)
"""

# high level stuff
include("Common.jl")

# low level containers
include("BoundaryConditions.jl")
include("PostProcessors.jl")
include("Sections.jl")
include("TimeSteppers.jl")
include("Domains.jl")

# solver
include("solvers/NonlinearSolvers.jl")

# include("Objective.jl")

include("Problems.jl")

# CLI
include("Main.jl")

end # module
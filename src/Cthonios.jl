module Cthonios

# TODO add exports

# dependencies
using ArgParse
using ConstitutiveModels
using DocStringExtensions
using Exodus
using FiniteElementContainers
using ForwardDiff
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
using SparseArrays
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

# Level 1 stuff
include("Common.jl")
include("Functions.jl")
include("Parsers.jl")
include("PostProcessors.jl")

# Level 2 stuff
include("boundary_conditions/BoundaryConditions.jl")
include("sections/Sections.jl")
include("TimeSteppers.jl")

# low level containers
# include("Backends.jl")
# include("BoundaryConditions.jl")
# include("PostProcessors.jl")
# include("Sections.jl")
# include("TimeSteppers.jl")
# include("Domains.jl")

include("domains/Domains.jl")

# solver
include("solvers/NonlinearSolvers.jl")


include("Problems.jl")

# Parsing last thing before main files


# CLI
include("Main.jl")

end # module
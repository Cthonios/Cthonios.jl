module Cthonios

# TODO add exports

# dependencies
using ArgParse
using Atomix
using ComponentArrays
using ConstitutiveModels
using DocStringExtensions
using Exodus
using FiniteElementContainers
using FunctionWrappers
using IterativeSolvers
using KernelAbstractions
using LinearAlgebra
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

# Level 3 stuff
include("domains/Domains.jl")
include("Mechanics.jl")
include("solvers/NonlinearSolvers.jl")

# Finally problems
include("Problems.jl")

# CLI
include("Main.jl")

end # module
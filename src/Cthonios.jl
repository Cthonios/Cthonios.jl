module Cthonios

# edpxort
export CthoniosCommon
export Domain
export EssentialBC

# export create_field

# dev stuff

# export TotalLagrangeSection


# dependencies
# using AbstractDifferentiation
using ConstitutiveModels
using DocStringExtensions
using Enzyme
using Exodus
using FiniteElementContainers
using FunctionWrappers
using IterativeSolvers
using LinearAlgebra
using LinearMaps
using Logging
using LoggingExtras
using Parameters
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

# CLI
include("Main.jl")

end # module
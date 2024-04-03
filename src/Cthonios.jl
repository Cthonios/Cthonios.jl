module Cthonios

# TODO add exports

# dependencies
using ArgParse
using ComponentArrays
using ConstitutiveModels
using DocStringExtensions
using Enzyme
using Exodus
using FiniteElementContainers
using FunctionWrappers
using IterativeSolvers
using KernelAbstractions
using LDLFactorizations
using LinearAlgebra
using LinearOperators
using Logging
using LoggingExtras
using Parameters
using Pkg
using Printf
using ReferenceFiniteElements
using SparseArrays
using SparseDiffTools
using StaticArrays
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
# include("Backends.jl")
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
include("mechanics/Mechanics.jl")
include("solvers/WarmStart.jl")
include("solvers/NonlinearSolvers.jl")

# Finally problems
include("problems/Problems.jl")

# CLI
include("Main.jl")

# Precompile the features you want to be fast in executables
# or just to minimize time to first action in the REPL
# @setup_workload begin
#   for dir in readdir("precompile")
#     cthonios_main(joinpath("precompile", dir, "precompile.yaml"), false, "CPU")
#   end
# end

end # module
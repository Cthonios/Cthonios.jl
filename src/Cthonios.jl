module Cthonios

# Re-exports
export DirichletBC
export InitialCondition
export PlaneStrain
export ThreeDimensional
export TimerOutput
export TimeStepper
export UnstructuredMesh
export energy
export residual
export stiffness
export @SVector

# Cthonios exports
export ContactPair

# objectives
export EigenObjective
export ExplicitDynamicsObjective
export ImplicitDynamicsObjective
export QuasiStaticObjective

# physics
export SolidMechanics

# qois
export QOIExtractor

# sims
export SingleDomainSimulation

# solver stuff
# export AMGPreconditioner
export CholeskyPreconditioner
export EigenSolver
export ExplicitSolver
export ImplicitSolver
export JacobiPreconditioner
export LLDLPreconditioner
export LUPreconditioner
export NoPreconditioner
export NoPredictor
export TangentPredictor
export TrustRegionSolver
export create_unknowns
export evolve!

import FiniteElementContainers as FEC
import FiniteElementContainers: AbstractField, BCBookKeeping
import KernelAbstractions as KA
import KernelAbstractions: CPU
# using AlgebraicMultigrid
using Arpack
using ConstitutiveModels
using DocStringExtensions
using FiniteElementContainers
using ForwardDiff
using Krylov
using LimitedLDLFactorizations
using LinearAlgebra
using NLopt
using Printf
using ReferenceFiniteElements
using SparseArrays
using StaticArrays
using StructArrays
using Tensors
using TimerOutputs


mutable struct DeveloperOptions
    use_condensed::Bool
    use_inplace_methods::Bool

    function DeveloperOptions()
        new(false, true)
    end
end

const DEV_OPTIONS = DeveloperOptions()

function set_developer_option!(key::Symbol, val)
    setfield!(DEV_OPTIONS, key, val)
    return nothing
end

# objectives
include("objectives/Objectives.jl")

# physics
# include("physics/Physics.jl")
include("Physics.jl")

# contact
include("contact/Contact.jl")

# simulations
# include("PostProcessor.jl")
include("simulations/Simulations.jl")

# solvers
include("solvers/Solvers.jl")

#
include("qoi_extractors/QOIExtractors.jl")
# include("sensitivities/Sensitivities.jl")

# optimizations
include("optimizations/Optimizations.jl")

end # module

module Cthonios

import FiniteElementContainers: AbstractField, BCBookKeeping
import KernelAbstractions as KA
import KernelAbstractions: CPU
using Arpack
using ConstitutiveModels
using DocStringExtensions
using FiniteElementContainers
using ForwardDiff
using Krylov
using LinearAlgebra
using NLopt
using Printf
using ReferenceFiniteElements
using SparseArrays
using StaticArrays
using StructArrays
using Tensors
using TimerOutputs
# using YAML

# RuntimeGeneratedFunctions.init(@__MODULE__)

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
# export NewtonSolver # name conflict with FiniteElementContainers
export ExplicitDynamicsObjective
export ImplicitDynamicsObjective
export QOIExtractor
export QuasiStaticObjective
export SingleDomainSimulation
export SolidMechanics
export TrustRegionSolver
export create_unknowns
export evolve!

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
include("PostProcessor.jl")
include("simulations/Simulations.jl")

# solvers
include("solvers/Solvers.jl")

#
include("qoi_extractors/QOIExtractors.jl")
# include("sensitivities/Sensitivities.jl")

# optimizations
include("optimizations/Optimizations.jl")

# include("cli/CLI.jl")

# methods defined in extensions
# function cthonios_main end

# function @main(ARGS::Vector{String})
#     cthonios_main(ARGS)
#     return 0
# end

end # module

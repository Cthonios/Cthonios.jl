module Cthonios

import FiniteElementContainers: AbstractField, BCBookKeeping
import KernelAbstractions as KA
import KernelAbstractions: CPU
using ArgParse
using Arpack
using ConstitutiveModels
using DocStringExtensions
using Enzyme
using Exodus
using FiniteElementContainers
using ForwardDiff
using Krylov
using LinearAlgebra
using NLopt
using Printf
using ReferenceFiniteElements
using RuntimeGeneratedFunctions
using SparseArrays
using StaticArrays
using StructArrays
using Tensors
using TimerOutputs
using YAML

RuntimeGeneratedFunctions.init(@__MODULE__)

# Re-exports
export DirichletBC
export PlaneStrain
export QuasiStaticObjective
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
export SingleDomainSimulation
export SolidMechanics
export TrustRegionSolver
export create_unknowns
export evolve!
export parameters

# objectives
include("objectives/Objectives.jl")

# physics
include("physics/Physics.jl")

# contact
include("contact/Contact.jl")

# simulations
include("PostProcessor.jl")
include("simulations/Simulations.jl")

# solvers
include("solvers/Solvers.jl")

#
include("qoi_extractors/QOIExtractors.jl")

# optimizations
include("optimizations/Optimizations.jl")

include("cli/CLI.jl")

# methods defined in extensions
function cthonios_main end

end # module

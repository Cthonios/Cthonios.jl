module Cthonios

using ConstitutiveModels
using DocStringExtensions
using Enzyme
using Exodus
using FiniteElementContainers
using IncompleteLU
using KernelAbstractions
using Krylov
using LinearAlgebra
using NLopt
using Printf
using RuntimeGeneratedFunctions
using SparseArrays
using StaticArrays
using TimerOutputs

RuntimeGeneratedFunctions.init(@__MODULE__)

# Re-exports
export DirichletBC
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
export QuadratureLevelObjective
export SingleDomainSimulation
export SolidMechanics
export TrustRegionSolver
export UnconstrainedObjective
export create_unknowns
export evolve!
export parameters

# utilities
include("Utils.jl")

# objectives
include("objectives/Objectives.jl")

# physics
include("physics/Physics.jl")

# simulations
include("simulations/Simulations.jl")

# solvers
include("solvers/Solvers.jl")

#
include("qoi_extractors/QOIExtractors.jl")

# optimizations
include("optimizations/Optimizations.jl")

# CLI
# include("Main.jl")
function cthonios_main end

end # module

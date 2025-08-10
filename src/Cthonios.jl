module Cthonios

using ArgParse
using ConstitutiveModels
using DocStringExtensions
using Enzyme
using Exodus
using FiniteElementContainers
using LinearAlgebra
using NLopt
using Printf
using StaticArrays
using TimerOutputs
using YAML

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
include("Main.jl")

end # module

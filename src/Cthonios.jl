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

include("contact/Contact.jl")
include("objectives/Objectives.jl")
include("Physics.jl")
include("simulations/Simulations.jl")
include("solvers/Solvers.jl")

end # module

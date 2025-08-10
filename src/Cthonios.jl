module Cthonios

using ArgParse
using ConstitutiveModels
using DocStringExtensions
# using EngineeringSketchPadWrapper
using Enzyme
using Exodus
using FiniteElementContainers
using LinearAlgebra
# using NLopt
using Optimization
using OptimizationNLopt
using Printf
using StaticArrays
using TimerOutputs
using YAML

# Re-exports
export FiniteElementContainers
export StaticArrays

export SingleDomainSimulation
export UnconstrainedObjective
export parameters

# small components that make up a more complex problem
# include("bcs/BoundaryConditions.jl")
# include("physics/Physics.jl")
# include("PostProcessors.jl")
# include("Sections.jl")

# # container for small components
# include("Domains.jl")

# # extra optional stuff related to constraints
# include("contact/Contact.jl")

# # time integration facilities
# include("integrators/Integrators.jl")

# # integrals over objective function kernels
# include("integrals/Integrals.jl")

# # different objectives
# include("objectives/Objectives.jl")

# # integrals over objective function kernels
# # include("integrals/Integrals.jl")

# # solvers
# include("solvers/Solvers.jl")

# # different problem types
# include("simulations/Simulations.jl")
# # include("SchwarzProblem.jl")

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

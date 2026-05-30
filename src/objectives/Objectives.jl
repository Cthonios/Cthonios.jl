"""
User facing
"""
abstract type AbstractObjective{F1} end
abstract type AbstractSolutionObjective{F1} <: AbstractObjective{F1} end

# new cache implementation
abstract type AbstractObjectiveCache{
    A, # Assembler type
    O  <: AbstractObjective,
    RT <: Number,
    RV <: AbstractArray{RT, 1}
} end

function FiniteElementContainers.create_field(o::AbstractObjectiveCache)
    return FiniteElementContainers.create_field(o.assembler)
end

function FiniteElementContainers.create_unknowns(o::AbstractObjectiveCache)
    return FiniteElementContainers.create_unknowns(o.assembler)
end

function assembler(o::AbstractObjectiveCache)
    return o.assembler
end

abstract type AbstractSolutionObjectiveCache{A, O, RT, RV} <: AbstractObjectiveCache{A, O, RT, RV} end

# function parameters(o::AbstractObjectiveCache)
#     return o.parameters
# end
# include("ContactObjective.jl")
# include("ConstrainedObjective.jl")
# include("DesignObjective.jl")
include("ExplicitDynamicsObjective.jl")
include("ImplicitDynamicsObjective.jl")
include("QuasiStaticObjective.jl")

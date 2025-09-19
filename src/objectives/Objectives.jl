"""
User facing
"""
abstract type AbstractObjective{F1 <: Function} end

abstract type AbstractObjectiveCache{
    O <: AbstractObjective,
    S, # Type me!
    T <: TimerOutput
} end

function FiniteElementContainers.create_field(o::AbstractObjectiveCache)
    return FiniteElementContainers.create_field(o.sim_cache.assembler)
end

function FiniteElementContainers.create_unknowns(o::AbstractObjectiveCache)
    return FiniteElementContainers.create_unknowns(o.sim_cache.assembler)
end

function assembler(o::AbstractObjectiveCache)
    return o.sim_cache.assembler
end

function parameters(o::AbstractObjectiveCache)
    return o.sim_cache.parameters
end

# new cache implementation
abstract type AbstractObjectiveCache2{
    A, # Assembler type
    O  <: AbstractObjective,
    P  <: FiniteElementContainers.Parameters,
    RT <: Number,
    RV <: AbstractArray{RT, 1}
} end

function FiniteElementContainers.create_field(o::AbstractObjectiveCache2)
    return FiniteElementContainers.create_field(o.assembler)
end

function FiniteElementContainers.create_unknowns(o::AbstractObjectiveCache2)
    return FiniteElementContainers.create_unknowns(o.assembler)
end

function assembler(o::AbstractObjectiveCache2)
    return o.assembler
end

function parameters(o::AbstractObjectiveCache2)
    return o.parameters
end

include("QuadratureLevelObjective.jl")

include("ImplicitDynamicsObjective.jl")
include("ImplicitDynamicsObjectiveNew.jl")
include("QuasiStaticObjective.jl")
include("QuasiObjectiveNew.jl")


include("UnconstrainedObjective.jl")

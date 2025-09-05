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

include("QuadratureLevelObjective.jl")
include("UnconstrainedObjective.jl")

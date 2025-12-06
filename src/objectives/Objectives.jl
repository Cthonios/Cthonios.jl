"""
User facing
"""
abstract type AbstractObjective{
    F1 <: Function
} end

# new cache implementation
abstract type AbstractObjectiveCache{
    A, # Assembler type
    O  <: AbstractObjective,
    # P  <: FiniteElementContainers.Parameters,
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

# function parameters(o::AbstractObjectiveCache)
#     return o.parameters
# end

# include("DesignObjective.jl")
# include("ImplicitDynamicsObjective.jl")
include("QuasiStaticObjective.jl")

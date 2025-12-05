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
    P  <: FiniteElementContainers.Parameters,
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

function parameters(o::AbstractObjectiveCache)
    return o.parameters
end

# abstract fall backs for quicker prototyping
# could eventually make more efficient
function gradient(cache::AbstractObjectiveCache, U, p, ::Val{:enzyme})
    dU = make_zero(U)
    dp = make_zero(p)
    dcache = make_zero(cache)
    autodiff(
        Reverse,
        value,
        Duplicated(cache, dcache),
        Duplicated(U, dU),
        Duplicated(p, dp)
    )
    dU, dp
end

# # stuff for adjoints
# function dresidual_du(o::AbstractObjective, U, p)
#     # just make copies for now, optimize later
#     dU = make_zero(U)
#     dp = make_zero(p)
    
# end


include("ImplicitDynamicsObjective.jl")
include("QuasiStaticObjective.jl")

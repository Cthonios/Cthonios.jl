"""
$(TYPEDEF)
```NF``` - Number of fields in this physics\n
```NP``` - Number of properties in this physics\n
```NS``` - Number of states in this physics
"""
abstract type AbstractPhysics{NF, NP, NS} end

"""
$(TYPEDSIGNATURES)
"""
num_fields(::AbstractPhysics{NF, NP, NS}) where {NF, NP, NS} = NF
"""
$(TYPEDSIGNATURES)
"""
num_properties(::AbstractPhysics{NF, NP, NS}) where {NF, NP, NS} = NP
"""
$(TYPEDSIGNATURES)
"""
num_states(::AbstractPhysics{NF, NP, NS}) where {NF, NP, NS} = NS

# implementations
include("Poisson.jl")
include("SolidMechanics.jl")

# exports
export Poisson

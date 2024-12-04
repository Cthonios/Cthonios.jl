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

# implementations as small kernels which operate
# at the quadratue level
include("kernels/Dynamics.jl")
include("kernels/Laplacian.jl")
include("kernels/StressDivergence.jl")
include("kernels/Source.jl")

# useful combinations of kernels which also operate
# at the quadrature level
include("Poisson.jl")
include("SolidDynamics.jl")
include("SolidMechanics.jl")

# wrappers that operate at the element level
include("ElementLevel.jl")

# exports
export Poisson
export SolidDynamics
export SolidMechanics

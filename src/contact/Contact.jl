# Organize contact as follows
# 1. enforcement strategy, e.g. node-face, face-face, mortar
# 2. contact formulation, e.g. penalty, AL
#
abstract type AbstractContactEnforcement end
abstract type AbstractContactFormulation end

function assemble_contact_matrix!(
    enforcement::AbstractContactEnforcement,
    formulation::AbstractContactFormulation,
    caches::NamedTuple,
    X, U
)
    for cache in values(caches)
        assemble_contact_matrix!(enforcement, formulation, cache, X, U)
    end
end

function assemble_contact_scalar!(
    energy,
    enforcement::AbstractContactEnforcement,
    formulation::AbstractContactFormulation,
    caches::NamedTuple,
    # caches,
    X, U
)
    for cache in values(caches)
        assemble_contact_scalar!(energy, enforcement, formulation, cache, X, U)
    end
end

function assemble_contact_vector!(
    vector,
    enforcement::AbstractContactEnforcement,
    formulation::AbstractContactFormulation,
    caches::NamedTuple,
    X, U
)
    for cache in values(caches)
        assemble_contact_vector!(vector, enforcement, formulation, cache, X, U)
    end
end

# only for 2D
# assumes....
# closest point projection distance
include("AD.jl")
include("ContactSurface.jl")
include("ContactPair.jl")

include("ClosestPointProjection.jl")
include("Search.jl")

include("PenaltyContact.jl")
include("Mortar.jl")


# TODO contact_stiffness_side_a_side_b
# need to be careful when facet element types differ

# used to construct these guys
function create_contact_pair_caches(
    mesh, dof, contact_pairs;
    kwargs...
)
    @assert num_dimensions(mesh) == 2 "Contact is only supported in 2D currently"
    syms = map(x -> Symbol("contact_pair_$(x.side_a)_$(x.side_b)"), contact_pairs)
    caches = map(x -> ContactPairCache(mesh, dof, :displ, x; kwargs...), contact_pairs)
    return NamedTuple{tuple(syms...)}(tuple(caches...))
end

# Organize contact as follows
# 1. enforcement strategy, e.g. node-face, face-face, mortar
# 2. contact formulation, e.g. penalty, AL
#
abstract type AbstractContactEnforcement end
abstract type AbstractContactFormulation end

# only for 2D
# assumes....
# closest point projection distance
include("ContactSurface.jl")

include("ContactPair.jl")

include("PenaltyContact.jl")

include("Mortar.jl")
include("Search.jl")

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
    enforcement::AbstractContactEnforcement,
    formulation::AbstractContactFormulation,
    caches::NamedTuple,
    X, U
)
    for cache in values(caches)
        assemble_contact_scalar!(enforcement, formulation, cache, X, U)
    end
end

function assemble_contact_vector!(
    enforcement::AbstractContactEnforcement,
    formulation::AbstractContactFormulation,
    caches::NamedTuple,
    X, U
)
    for cache in values(caches)
        assemble_contact_vector!(enforcement, formulation, cache, X, U)
    end
end

function contact_residual_side_a(
    ::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.gradient(
        z -> contact_energy(
            X_a, X_b, z, U_b, 
            max_overlap_dist, rel_smoothing_size
        ), U_a
    )
end

function contact_residual_side_b(
    ::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.gradient(
        z -> contact_energy(
            X_a, X_b, U_a, z, 
            max_overlap_dist, rel_smoothing_size
        ), U_b
    )
end

function contact_stiffness_side_a_side_a(
    ::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.hessian(
        z -> contact_energy(
            X_a, X_b, z, U_b, 
            max_overlap_dist, rel_smoothing_size
        ), U_a
    )
end

function contact_stiffness_side_b_side_b(
    ::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.hessian(
        z -> contact_energy(
            X_a, X_b, U_a, z, 
            max_overlap_dist, rel_smoothing_size
        ), U_b
    )
end

function contact_stiffness_side_a_side_b(
    ::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.jacobian(
        z -> ForwardDiff.gradient(
            contact_residual_side_a(
                X_a, X_b, U_a, z, 
                max_overlap_dist, rel_smoothing_size
            ), U_b
        )
    )
end


# TODO contact_stiffness_side_a_side_b
# need to be careful when facet element types differ

# used to construct these guys
function create_contact_pair_caches(mesh, contact_pairs)
    @assert num_dimensions(mesh) == 2 "Contact is only supported in 2D currently"
    syms = map(x -> Symbol("contact_pair_$(x.side_a)_$(x.side_b)"), contact_pairs)
    caches = map(x -> ContactPairCache(mesh, x), contact_pairs)
    return NamedTuple{tuple(syms...)}(tuple(caches...))
end

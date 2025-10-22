struct ContactPair
    side_a::Symbol
    side_b::Symbol
end

function ContactPair(side_a::String, side_b::String)
    return ContactPair(Symbol(side_a), Symbol(side_b))
end

struct ContactPairCache{
    IT <: Integer,
    # RT <: Number,
    IM <: AbstractMatrix{IT},
    IV <: AbstractVector{IT},
    # RV <: AbstractVector{RT},
    DriverSurf <: ContactSurface{IT, IM, IV, <:ReferenceFE},
    FollowSurf <: ContactSurface{IT, IM, IV, <:ReferenceFE},
    # MaxNeighbors, NNPF, NQ
    MaxNeighbors
}
    side_a::DriverSurf
    side_b::FollowSurf
    #
    interactions::H1Field{IT, IV, MaxNeighbors}
    #
    # closest_facets::L2QuadratureField{IT, IV, NNPF, NQ}
    # facet_weights::L2QuadratureField{RT, RV, 1, NQ}
    #
    self_contact::Bool
end

function ContactPairCache(
    mesh, contact_pair::ContactPair;
    max_neighbors=1,
    q_order=1,
    self_contact=false
)
    side_a = ContactSurface(mesh, contact_pair.side_a, q_order)
    side_b = ContactSurface(mesh, contact_pair.side_b, q_order)
    return ContactPairCache(
        side_a, side_b;
        max_neighbors=max_neighbors,
        q_order=q_order,
        self_contact=self_contact
    )
end


function ContactPairCache(
    side_a::ContactSurface, side_b::ContactSurface;
    max_neighbors=4,
    q_order=1,
    self_contact=false
)
    # setup with interactions not initialized but sized correctly
    # n_driver_side_nodes = size(side_a.side_nodes, 2)
    n_sides = size(side_b.side_nodes, 2)
    interactions = H1Field(zeros(Int, max_neighbors, n_sides))

    # NNPF = num_vertices(surface_element(side_a.ref_fe.element))
    # closest_faces = L2QuadratureField(zeros(Int, NNPF, num_quadrature_points(side_a.ref_fe), n_sides))
    # facet_weights = L2QuadratureField(zeros(1, num_quadrature_points(side_a.ref_fe), n_sides))
    # closest_faces = nothing
    # facet_weights = nothing

    # return ContactPairCache(side_a, side_b, interactions, closest_faces, facet_weights, self_contact)
    return ContactPairCache(side_a, side_b, interactions, self_contact)
end

function Base.show(io::IO, cache::ContactPairCache)
    println(io, "ContactPairCache:")
    show(io, cache.side_a; tab="  ")
    show(io, cache.side_b; tab="  ")
end

function get_interaction(cache::ContactPairCache, n::Int, facet_id::Int)
    return cache.interactions[n, facet_id]
end

function get_interaction_facet(cache::ContactPairCache, n::Int, facet_id::Int)
    NN = num_vertices(surface_element(cache.side_a.ref_fe.element))
    side_a = cache.side_a.side_nodes
    return SVector{NN, Int}(@views side_a[:, get_interaction(cache, n, facet_id)])
end

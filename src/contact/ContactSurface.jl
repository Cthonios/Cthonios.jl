struct ContactSurface{
    IT <: Integer,
    IM <: AbstractMatrix{IT},
    IV <: AbstractVector{IT},
    RefFe <: ReferenceFE
}
    ref_fe::RefFe
    side_nodes::IM
    sides::IV
end

function ContactSurface(mesh, dof, var_name, sset_name::Symbol, q_order::Int)
    bk = BCBookKeeping(mesh, dof, var_name; sset_name=sset_name)

    # ensure there's only one block in each surface
    # and that the element types are equivalent
    @assert length(unique(bk.blocks)) == 1

    # grab the block to set up the reference FE
    block = read_set(mesh.mesh_obj.mesh_obj, Block, bk.blocks[1])
    ref_fe = ReferenceFE(block, q_order)

    # reshape side nodes
    nnps = num_vertices(surface_element(ref_fe.element))
    side_nodes = reshape(bk.side_nodes, nnps, length(bk.sides))
    sides = bk.sides
    return ContactSurface(ref_fe, side_nodes, sides)
end

function Base.show(io::IO, surf::ContactSurface; tab="")
    println(io, "$(tab)ContactSurface:")
    println(io, "$(tab)  Number of facets     = $(length(surf.sides))")
    println(io, "$(tab)  Surface element type = $(surface_element(surf.ref_fe.element))")
end

# returns a matrix with all the nodal coordinates in the facet
@inline function facet_coordinates(surf::ContactSurface, X, U, id::Int)
    RT = eltype(X)
    ND = size(X, 1)
    NNPF = 2 # hardcoded for now TODO

    coords = SMatrix{ND, NNPF, RT, ND * NNPF}(
        @views X[:, surf.side_nodes[:, id]]
    )
    disps = SMatrix{ND, NNPF, RT, ND * NNPF}(
        @views U[:, surf.side_nodes[:, id]]
    )
    return coords + disps
end

# returns a vector for a particular node's coordinates in the facet
@inline function facet_coordinates(surf::ContactSurface, X, U, node::Int, facet::Int)
    RT = eltype(X)
    ND = size(X, 1)
    NNPF = 2 # hardcoded for now TODO

    coords = SVector{ND, RT}(
        @views X[:, surf.side_nodes[node, facet]]
    )
    disps = SVector{ND, RT}(
        @views U[:, surf.side_nodes[node, facet]]
    )
    return coords + disps
end

# returns flat things
@inline function facet_field(surf::ContactSurface, X, id::Int)
    RT = eltype(X)
    ND = size(X, 1)
    NNPF = 2 # hardcoded for now

    # coords = SMatrix{ND, NNPF, RT, ND * NNPF}(
    #     @views X[:, surf.side_nodes[:, id]]
    # )
    coords = SVector{ND * NNPF, RT}(
        @views X[:, surf.side_nodes[:, id]]
    )
    return coords
end

@inline function facet_field(surf::ContactSurface, X, node::Int, facet::Int)
    RT = eltype(X)
    ND = size(X, 1)
    NNPF = 2 # hardcoded for now

    coords = SVector{ND, RT}(
        @views X[:, surf.side_nodes[node, facet]]
    )
    return coords
end

# below X_f is an already calculated facet coordinate
# using above
@inline function facet_quadrature_point_coordinates(surf::ContactSurface, X_f, q::Int, side::Int)
    interps = MappedSurfaceInterpolants(surf.ref_fe, X_f, q, side)
    return interps.X_q
end

# only for edge
# TODO dispatch on surface element
# currently specialized to linear edge
# case for when flat vecs come in, better for AD
@inline function facet_normal(X_f::SVector{4, RT}) where RT <: Number
    X_f = SMatrix{2, 2, RT, 4}(X_f.data)
    return facet_normal(X_f)
end

@inline function facet_normal(X_f::SMatrix{2, 2, RT, 4}) where RT <: Number
    tangent = SVector{2, RT}(
        X_f[1, 2] - X_f[1, 1],
        X_f[2, 2] - X_f[2, 1]
    )
    normal = SVector{2, RT}(tangent[2], -tangent[1])
    return normal / norm(normal)
end

@inline function facet_normal_from_a(X_a, X_b)
    return facet_normal(X_a)
end

@inline function facet_normal_from_b(X_a, X_b)
    return -facet_normal(X_b)
end

function _compute_gap(edge_a, xi_a, edge_b, n)
    # x_a = edge_a[SVector{2, Int}(1, 2)] * (1. - xi_a) + 
    #       edge_a[SVector{2, Int}(3, 4)] * xi_a
    x_a = edge_a[SVector{2, Int}(1, 2)] * (1. - xi_a) + 
          edge_a[SVector{2, Int}(3, 4)] * xi_a
    _, g_b = _compute_xi(x_a, edge_b, n)
    return g_b
end

function _compute_xi(X_a, X_b, n)
    # solve xa - xb(xi) + g * normal = 
    # RT = eltype(X_b)
    RT = promote_type(eltype(X_a), eltype(X_b), eltype(n))
    # M = SMatrix{2, 2, RT, 4}(X_b[2], n[1], X_b[1], n[2])
    # M = SMatrix{2, 2, RT, 4}(X_b[2], X_b[1], n[1], n[2])
    # M = SMatrix{2,2,RT}(
    #     (X_b[1] - X_b[2]), 
    #     # n[1],
    #     (X_b[3] - X_b[4]), 
    #     # n[2]
    #     n[1], n[2]
    # )
    M = SMatrix{2, 2, RT, 4}(
        (X_b[1] - X_b[3]), 
        # n[1],
        (X_b[2] - X_b[4]), 
        # (X_b[SVector{2, Int}(1, 2)] - X_b[SVector{2, Int}(3, 4)])...,
        n[1], 
        n[2]
    )
    # display(M)
    detM = det(M)
    if abs(detM) < 1.e-13 || dot(n, n) < 0.9
        M = SMatrix{2, 2, RT, 4}(1., 0., 0., 1.)
    end

    r = X_b[SVector{2, Int}(1, 2)] - X_a
    # display(r)
    xig = M \ r
    return xig[1], xig[2]
end

function _smooth_linear(xi, l)
    if xi < l
        temp = 0.5 * xi * xi / l
    else
        if xi > 1. - l
            temp = 1.0 - l - 0.5 * (1.0 - xi) * (1.0 - xi) / l
        else
            temp = xi - 0.5 * l
        end
    end

    return temp / (1. - l)
end

# TODO need to specialize on element type
function intersection(
    X_a, X_b, n, rel_smoothing_size
)
    # display(X_a)
    # display(X_b)
    # display(n)
    RT = promote_type(eltype(X_a), eltype(X_b), eltype(n))
    xi_g_b_right = _compute_xi(X_a[SVector{2, Int}(1, 2)], X_b, n)
    xi_g_a_right = _compute_xi(X_b[SVector{2, Int}(3, 4)], X_a, -n)
    xi_g_b_left = _compute_xi(X_a[SVector{2, Int}(3, 4)], X_b, n)
    xi_g_a_left = _compute_xi(X_b[SVector{2, Int}(1, 2)], X_a, -n)

    xi_g_b_right_inside = xi_g_b_right[1] <= 1.
    xi_g_b_left_inside = xi_g_b_left[1] >= 0.

    if xi_g_b_right_inside
        temp_1 = 0.0
    else
        temp_1 = xi_g_a_right[1]
    end

    if xi_g_b_left_inside
        temp_2 = 1.
    else
        temp_2 = xi_g_a_left[1]
    end

    xi_minmax_in_a = SVector{2, RT}(temp_1, temp_2)
    xi_minmax_in_a = map(x -> max(0., min(1., x)), xi_minmax_in_a)
    xi_minmax_in_a = _smooth_linear.(xi_minmax_in_a, (rel_smoothing_size,))

    gaps_a = map(x -> _compute_gap(X_a, x, X_b, n), xi_minmax_in_a)

    return xi_minmax_in_a, gaps_a
end

# TODO incorporate self contact into search or
# dispatch on it somehow

# only works for linear edge elements
# function _edge_is_self(cache, a, b)
#     side_a = cache.side_a.side_nodes
#     side_b = cache.side_b.side_nodes
#     return side_a[1, a] == side_b[1, b] && \
#         side_b[2, a] == side_b[2, b]
# end

# # only works for linear edge elements
# function _edges_are_adjacent_non_pacman(cache, X, a, b)
#     side_a = cache.side_a.side_nodes
#     side_b = cache.side_b.side_nodes
#     adj_case_1 = side_a[2, a] == side_b[1, b]
#     adj_case_2 = side_a[1, a] == side_b[2, b]


# end

@inline function _min_distance_between_edges(cache, X, U, edge_1_index, edge_2_index)
    min_distance = 1e12

    for dn in axes(cache.side_a.side_nodes, 1)
        for fn in axes(cache.side_b.side_nodes, 1)
            u_d = facet_coordinates(cache.side_a, X, U, dn, edge_1_index)
            u_f = facet_coordinates(cache.side_b, X, U, fn, edge_2_index)
            temp_dist = dot(u_d - u_f, u_d - u_f)
            min_distance = min(min_distance, temp_dist)
        end
    end
    return min_distance
end

function update_nearest_neighbors!(
    cache::ContactPairCache,
    X, U # eventually make distance func an input
)
    @assert size(X) == size(U)

    distances = Vector{eltype(X)}(undef, 0)
    closest_ids = Vector{Int}(undef, 0)

    max_neighbors = size(cache.interactions, 1)

    for b in axes(cache.side_b.side_nodes, 2)
        resize!(distances, 0)
        for a in axes(cache.side_a.side_nodes, 2)
            if cache.self_contact
                @assert false "TODO implement self contact corrections"
            else
                min_dist = _min_distance_between_edges(cache, X, U, a, b)
            end
            push!(distances, min_dist)
        end
        # @show distances
        resize!(closest_ids, length(distances))
        sortperm!(closest_ids, distances)
        sort!(distances)

        for n in axes(closest_ids, 1)
            if !isfinite(distances[n])
                closest_ids[n] = -1
            end
        end

        @views cache.interactions[:, b] = closest_ids[1:max_neighbors]
    end

    return nothing
end

function update_nearest_neighbors!(cache::NamedTuple, X, U)
    for pair in values(cache)
        update_nearest_neighbors!(pair, X, U)
    end
    return nothing
end

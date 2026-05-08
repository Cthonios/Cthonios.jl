function _compute_normal(edge_coords)
    RT = eltype(edge_coords)
    tangent = SVector{2, RT}(
        edge_coords[1, 2] - edge_coords[1, 1],
        edge_coords[2, 2] - edge_coords[2, 1]
    )
    normal = SVector{2, eltype(edge_coords)}(tangent[2], -tangent[1])
    return normal / norm(normal)
end

function _cpp_line(edge, p)
    # @show edge
    a = edge[SVector{2, Int}(1, 2)]
    b = edge[SVector{2, Int}(3, 4)]
    v = b - a
    t = -dot(v, a - p) / dot(v, v)

    if t < 0.
        t = 0
    elseif t > 1.
        t = 1.
    end

    return (1. - t) * a + t * b, t
end

function _cpp_distance(edge, p)
    normal = _compute_normal(edge)
    cpp_pt, t = _cpp_line(edge, p)
    d = dot(normal, p - cpp_pt)
    sgn = sign(d)

    if sgn == 0.0
        sgn = 1.
    elseif sgn == -0.0
        sgn = 1.
    end

    if t < 0.
        temp = edge[SVector{2, Int}(1, 2)]
        d = sgn * sqrt(dot(temp - p))
    elseif t > 1.
        temp = edge[SVector{2, Int}(3, 4)]
        d = sgn * sqrt(dot(temp - p))
    end

    return d
end

function _compute_closest_distance_to_each_side(
    contact_caches, X, U
)
    # loop over contact pairs
    cpps = Float64[]
    for cpair in values(contact_caches)
        # loop over interactions in current contact pair
        for n in axes(cpair.interactions, 2)
            side = cpair.side_b.sides[n]
            X_I = facet_coordinates(cpair.side_b, X, U, n)
            # loop over quadrature points
            cpp = 1e64 # some really big number
            for q in 1:num_quadrature_points(cpair.side_b.ref_fe)
                X_Q = facet_quadrature_point_coordinates(cpair.side_b, X_I, q, side)
                # loop over neighbors in interaction list
                for m in axes(cpair.interactions, 1)
                    X_M = facet_coordinates(cpair.side_a, X, U, cpair.interactions[m, n])
                    cpp = min(cpp, _cpp_distance(X_M, X_Q))
                end
                # display(X_Q)
            end
            # @show cpp
            push!(cpps, cpp)
        end
    end
    cpps
end

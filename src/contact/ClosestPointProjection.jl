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
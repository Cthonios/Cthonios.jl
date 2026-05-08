struct PenaltyContact <: AbstractContactFormulation
end

function contact_energy(
    ::PenaltyContact,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    x_a = X_a + U_a
    x_b = X_b + U_b

    normal = facet_normal_from_b(x_a, x_b)
    n_a = facet_normal(x_a)
    n_b = facet_normal(x_b)

    pos_check = -dot(n_a, n_b)
    scaling = ifelse(pos_check > zero(eltype(X_a)), pos_check, zero(eltype(X_a)))
    x_b = ifelse(pos_check < zero(eltype(X_a)), x_a, x_b)

    # display(x_a)
    xi_a, g_a = intersection(x_a, x_b, normal, rel_smoothing_size)
    g = scaling * g_a + max_overlap_dist * ones(typeof(g_a))

    # TODO move below to method
    l_a = norm(X_a[SVector{2, Int}(1, 2)] - X_a[SVector{2, Int}(3, 4)])
    integral = integrate_gap(xi_a, g, max_overlap_dist)
    return l_a * integral * _smooth_heaviside_at_zero(scaling, 1.)
end

function _smooth_heaviside_at_zero(x, eps)
    r = x / eps
    # return jnp.where(x < eps, -2*r*r*r+3*r*r, 1.0)
    if x < eps
        return -2. * r * r * r + 3. * r * r
    else
        return 1.
    end
end

@inline function integrate_normalized_gap(
    ξ::AbstractVector{T},
    g::AbstractVector{T}
) where {T}
    dξ = ξ[2] - ξ[1]
    g0_large = g[1] >= one(eltype(g))
    g1_large = g[2] >= one(eltype(g))
    dξ_invalid = (dξ <= 1e-14) || (g0_large && g1_large)
    
    a = (g[2] - g[1]) / ifelse(dξ_invalid, one(eltype(g)), dξ)
    g = ifelse(g0_large, typeof(g)(one(eltype(g)), g[2]), g)
    g = ifelse(g1_large, typeof(g)(g[1], one(eltype(g))), g)

    a_is_zero = abs(g[1] - g[2]) < 1e-14
    a = ifelse(a_is_zero, one(eltype(a)), a)
    g0_old = g[1]
    g1_old = g[2]

    ξ = ifelse(
        g0_large && !a_is_zero, 
        typeof(ξ)((1. - g1_old) / a + ξ[2], ξ[2]),
        ξ
    )
    ξ = ifelse(
        g1_large && !a_is_zero, 
        typeof(ξ)(ξ[1], (1. - g0_old) / a + ξ[1]),
        ξ
    )
    dξ = ξ[2] - ξ[1]

    ξ0 = ξ[1]
    ξ1 = ξ[2]

    safe_log = if g[2] / g[1] < 0.
        0.
    else
        log(g[2] / g[1])
    end

    intgl = safe_log / a +
            (0.5 * a * (ξ1 * ξ1 - ξ0 * ξ0)) +
            g[1] * ξ1 - g[2] * ξ0 -
            2. * dξ
    
    gbar = 0.5 * (g[1] + g[2])
    dg = g[2] - g[1]
    intgl = ifelse(
        abs(g[1] - g[2]) < 2.e-4 * gbar,
        (1.0 / gbar + gbar - 2.0 + dg * dg / (12.0 * gbar * gbar * gbar)) * dξ, 
        intgl
    )
    any_neg = any(x -> x <= 0, g)
    intgl = ifelse(any_neg, Inf, intgl)
    intgl = ifelse(dξ_invalid, 0. * intgl, intgl)
    return intgl
end

@inline function integrate_gap(
    xi::AbstractVector{T},
    g::AbstractVector{T},
    delta
) where {T}
    return integrate_normalized_gap(xi, g ./ delta)
end

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
    integral = _integrate_gap(xi_a, g, max_overlap_dist)
    return l_a * integral * _smooth_heaviside_at_zero(scaling, 1.)
end

function _integrate_gap(xi, g, delta)
    return _integrate_normalized_gap(xi, g / delta)
end

# need these otherwise we allocate for some weird reason
# maybe we should use @setindex tools
_set_first(g::SVector{2, RT}, v::RT) where RT = SVector{2, RT}(v, g[2])
_set_second(g::SVector{2, RT}, v::RT) where RT = SVector{2, RT}(g[1], v)

function _integrate_normalized_gap(xi::SVector{2, RT}, g::SVector{2, RT}) where RT <: Number
    g
    dxi = xi[2] - xi[1]
    g1_large = g[1] >= one(RT)
    g2_large = g[2] >= one(RT)
    dxi_invalid = (dxi <= 1e-14) || (g1_large && g2_large)

    # a = (g[2] - g[1]) / jnp.where(dxi_invalid, 1.0, dxi)
    a = (g[2] - g[1])
    if !dxi_invalid
        a = a / dxi
    end

    # g = jnp.where(g0_large, g.at[0].set(1.0), g)
    if g1_large
        # g = SVector{2, RT}(1., g[2])
        g = _set_first(g, one(eltype(g)))
    end
    # g = jnp.where(g1_large, g.at[1].set(1.0), g)
    if g2_large
        # g = SVector{2, RT}(g[1], 1.)
        g = _set_second(g, one(eltype(g)))
    end

    a_is_zero = abs(g[1] - g[2]) < 1e-14
    # a = jnp.where(a_is_zero, 1.0, a)
    if a_is_zero
        a = 1.
    end

    g1_old = g[1]
    g2_old = g[2]

    # xi = jnp.where(g0_large & ~a_is_zero,
    #                xi.at[0].set((1.-g1_old) / a + xi[1]),
    #                xi)

    if g1_large && !a_is_zero
        xi = typeof(xi)((1. - g2_old) / a + xi[2], xi[2])
    end

    # xi = jnp.where(g1_large & ~a_is_zero,
    #                xi.at[1].set((1.0-g0_old) / a + xi[0]),
    #                xi)
    if g2_large && !a_is_zero
        xi = typeof(xi)(xi[1], (1. - g1_old) / a + xi[1])
    end

    dxi = xi[2] - xi[1]

    xi1 = xi[1]
    xi2 = xi[2]

    intgl = log(g[2] / g[1]) / a + (0.5 * a * (xi2 * xi2 - xi1 * xi1) + g[1] * xi2 - g[2] * xi1) - 2.0 * dxi

    gbar = 0.5 * (g[1] + g[2])
    dg = g[2] - g[1]
    # intgl = jnp.where(jnp.abs(g[0]-g[1]) < 2.0e-4 * gbar, (1.0/gbar + gbar - 2.0 + dg * dg / (12.0 * gbar*gbar*gbar)) * dxi, intgl)

    if abs(g[1] - g[2]) < 2.e-4 * gbar
        intgl = (1.0 / gbar + gbar - 2.0 + dg * dg / (12.0 * gbar * gbar * gbar)) * dxi
    end

    # any_neg = any(g <= 0)
    any_neg = any(x -> x <= zero(eltype(g)), g)
    # intgl = jnp.where(any_neg, jnp.inf, intgl)
    if any_neg
        intgl = Inf
    end
    # intgl = jnp.where(dxi_invalid, 0.0 * intgl, intgl)
    if dxi_invalid
        intgl = 0. * intgl
    end

    return intgl
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

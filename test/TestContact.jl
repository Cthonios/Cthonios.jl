using Exodus
using ReferenceFiniteElements
using Test

function single_quad_mesh_helper(mn)
    coords = [
        0. 1. 1. 0.;
        0. 0. 1. 1.
    ]

    # define edges as pairs of node indices
    edges_a = [
        1 2;
        2 3
    ]

    edges_b = [
        3 4;
        4 1
    ]

    disp = zeros(2, 4)

    ref_fe = ReferenceFE(Edge2{Lagrange, 1}())

    side_a = Cthonios.ContactSurface(ref_fe, edges_a, ones(Int, 2))
    side_b = Cthonios.ContactSurface(ref_fe, edges_b, ones(Int, 2))

    pair = Cthonios.ContactPairCache(side_a, side_b; max_neighbors=mn)

    return coords, pair
end

function two_block_tri3_mesh_helper(mn)
    mesh_file = Base.source_dir() * "/two_block_tri3_mesh.exo"
    mesh = UnstructuredMesh(mesh_file)
    fspace = FunctionSpace(mesh, H1Field, Lagrange)
    func = VectorFunction(fspace, :displ)
    dof = DofManager(func; use_condensed=true)
    pair = Cthonios.ContactPair("interface_top", "interface_bottom")
    cache = Cthonios.ContactPairCache(mesh, dof, :displ, pair; max_neighbors=mn)
    return mesh.nodal_coords, cache
end

function test_contact_pair_constructor()
    pair = Cthonios.ContactPair("cs1", "cs2")
    @test pair.side_a == :cs1
    @test pair.side_b == :cs2
end

function test_contact_surface_constructor_and_show()
    mesh_file = Base.source_dir() * "/two_block_tri3_mesh.exo"
    mesh = UnstructuredMesh(mesh_file)
    fspace = FunctionSpace(mesh, H1Field, Lagrange)
    func = VectorFunction(fspace, :displ)
    dof = DofManager(func; use_condensed=true)
    surf = Cthonios.ContactSurface(mesh, dof, :displ, Symbol("interface_top"), 1)
    @show surf
    # TODO do some actual testing other than just constructors
end

function test_neareast_neighbor_search_single_quad()
    X, pair = single_quad_mesh_helper(1)
    U = zeros(size(X))
    Cthonios.update_nearest_neighbors!(pair, X, U)
    @test pair.interactions[1, 1] == 2
    @test pair.interactions[1, 2] == 1

    X, pair = single_quad_mesh_helper(2)
    Cthonios.update_nearest_neighbors!(pair, X, U)
    @test pair.interactions[1, 1] == 2
    @test pair.interactions[2, 1] == 1
    @test pair.interactions[1, 2] == 1
    @test pair.interactions[2, 2] == 2
end

function test_neareast_neighbor_facet_field_single_quad()
    X, pair = single_quad_mesh_helper(1)
    U = zeros(size(X))

    side_a = pair.side_a

    # test getting coordinates and normal on surface
    X_a = Cthonios.facet_field(side_a, X, 1)
    n_a = Cthonios.facet_normal(X_a)
    @test n_a == SVector{2, Float64}(0., -1.)

    X_a = Cthonios.facet_field(side_a, X, 2)
    n_a = Cthonios.facet_normal(X_a)
    @test n_a == SVector{2, Float64}(1., 0.)

    X_a = Cthonios.facet_field(side_a, X, 1, 1)
    @test X_a == SVector{2, Float64}(0., 0.)
    X_a = Cthonios.facet_field(side_a, X, 2, 1)
    @test X_a == SVector{2, Float64}(1., 0.)
    
    X_a = Cthonios.facet_field(side_a, X, 1, 2)
    @test X_a == SVector{2, Float64}(1., 0.)
    X_a = Cthonios.facet_field(side_a, X, 2, 2)
    @test X_a == SVector{2, Float64}(1., 1.)
    
    side_b = pair.side_b

    # test getting coordinates and normal on surface
    X_b = Cthonios.facet_field(side_b, X, 1)
    n_b = Cthonios.facet_normal(X_b)
    @test n_b == SVector{2, Float64}(0., 1.)

    X_b = Cthonios.facet_field(side_b, X, 2)
    n_b = Cthonios.facet_normal(X_b)
    @test n_b == SVector{2, Float64}(-1., 0.)

    X_b = Cthonios.facet_field(side_b, X, 1, 1)
    @test X_b == SVector{2, Float64}(1., 1.)
    X_b = Cthonios.facet_field(side_b, X, 2, 1)
    @test X_b == SVector{2, Float64}(0., 1.)
    
    X_b = Cthonios.facet_field(side_b, X, 1, 2)
    @test X_b == SVector{2, Float64}(0., 1.)
    X_b = Cthonios.facet_field(side_b, X, 2, 2)
    @test X_b == SVector{2, Float64}(0., 0.)
end

# the below test assumes that the the contact surfaces
# are identical in element size and count
function test_neareast_neighbor_search_two_block_tri3_mesh(mn)
    X, pair = two_block_tri3_mesh_helper(mn)
    U = zeros(size(X))
    Cthonios.update_nearest_neighbors!(pair, X, U)
    side_a = reinterpret(SVector{2, Int}, pair.side_a.side_nodes)
    side_b = reinterpret(SVector{2, Int}, pair.side_b.side_nodes)

    for n in axes(pair.side_b.side_nodes, 2)
        x_b = SMatrix{2, 2, Float64, 4}(sortslices(X[:, pair.side_b.side_nodes[:, n]], dims=2))
        x_as = [X[:, Cthonios.get_interaction_facet(pair, m, n)] for m in axes(pair.interactions, 1)]
        x_as = map(x -> sortslices(x, dims=2), x_as)
        x_as = map(x -> SMatrix{2, 2, Float64, 4}(x), x_as)
        found_one = false
        for x_a in x_as
            for a in eachcol(x_a)
                if any(isapprox(SVector{2, Float64}(a), col; atol=5e-4) for col in eachcol(x_b))
                    found_one = true
                    break
                end
            end
            @test found_one
        end
    end
end

function test_integral_in_overlap()
    penalty_length = 0.1
    edge_smoothing = 0.2
    ξ = SVector{2, Float64}(0.2, 0.9)
    δ = 0.5
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(-0.1, 0.1), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(0.1, -0.1), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(-0.1, -0.1), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(-0.1, -0.2), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(-0.1, -0.1 + 1e-12), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(-0.2, 0.6), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(0.7, -0.3), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(0.7, 0.0), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(0.0, 0.7), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(0.2, 0.0), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(0.0, 0.3), δ) ≈ Inf
    @test Cthonios.integrate_gap(ξ, SVector{2, Float64}(0.0, 0.0), δ) ≈ Inf
end

# function _integrate_gap_numeric(xi, g, delta)
#     N = 10_000
#     xig = range(0.5/N, stop = 1.0 - 0.5/N, length = N)

#     dxi = xi[2] - xi[1]
#     w = dxi / N

#     gap(x) = g[1] + x * (g[2] - g[1])

#     p(x) = begin
#         v = gap(x)
#         v < delta ? (v / delta + delta / v - 2) : 0.0
#     end

#     # return sum(p, xig) * w
#     return sum(map(p, xig)) * w
# end

function _integrate_gap_numeric(xi, g, delta)
    N = 10_000

    ξ0, ξ1 = xi
    g0, g1 = g

    # enforce positive orientation (important!)
    if ξ1 < ξ0
        ξ0, ξ1 = ξ1, ξ0
        g0, g1 = g1, g0
    end

    dξ = ξ1 - ξ0
    w  = dξ / N

    # midpoint rule on [0,1]
    xg = range(0.5/N, stop = 1.0 - 0.5/N, length = N)

    function gap_at_x(x)
        ξ = ξ0 + x * dξ
        return g0 + (ξ - ξ0) / dξ * (g1 - g0)
    end

    function p(x)
        v = gap_at_x(x)
        if v <= 0
            return Inf
        elseif v < delta
            return v / delta + delta / v - 2
        else
            return 0.0
        end
    end

    return sum(p(x) for x in xg) * w
end

function test_integral_both_in_contact()
    ξ = SVector{2, Float64}(0.2, 0.9)
    g1 = SVector{2, Float64}(0.1, 0.3)
    g2 = SVector{2, Float64}(0.3, 0.1)
    δ = 0.5
    integral_1 = Cthonios.integrate_gap(ξ, g1, δ)
    integral_2 = Cthonios.integrate_gap(ξ, g2, δ)
    @test integral_1 ≈ integral_2

    integral_3 = _integrate_gap_numeric(ξ, g1, δ)
    integral_4 = _integrate_gap_numeric(ξ, g2, δ)

    @test isapprox(integral_1, integral_3, rtol = 1e-7)
    @test isapprox(integral_2, integral_4, rtol = 1e-7)
end

function test_integral_both_equal_in_contact()
    ξ = SVector{2, Float64}(0.2, 0.9)
    g = SVector{2, Float64}(0.3, 0.3)
    δ = 0.5
    integral_1 = _integrate_gap_numeric(ξ, g, δ)
    integral_2 = Cthonios.integrate_gap(ξ, g, δ)
    @test integral_1 ≈ integral_2
    # g[2] = g[2] + 1e-8
    g = SVector{2, Float64}(0.3, 0.3 + 1e-8)
    integral_3 = Cthonios.integrate_gap(ξ, g, δ)
    @test isapprox(integral_2, integral_3, atol=1e-7)
end

function test_integral_both_out_of_contact()
    ξ = SVector{2, Float64}(0.2, 0.9)
    g1 = SVector{2, Float64}(0.6, 0.8)
    g2 = SVector{2, Float64}(0.8, 0.6)
    δ = 0.5
    integral_1 = Cthonios.integrate_gap(ξ, g1, δ)
    integral_2 = Cthonios.integrate_gap(ξ, g2, δ)
    @test integral_1 ≈ 0.0
    @test integral_2 ≈ 0.0
    integral_3 = _integrate_gap_numeric(ξ, g1, δ)
    @test isapprox(integral_1, integral_3, atol=1e-7)
end

function test_integral_both_equal_out_of_contact()
    ξ = SVector{2, Float64}(0.2, 0.9)
    g1 = SVector{2, Float64}(0.8, 0.8)
    g2 = SVector{2, Float64}(0.8, 0.8)
    δ = 0.5
    integral_1 = Cthonios.integrate_gap(ξ, g1, δ)
    integral_2 = Cthonios.integrate_gap(ξ, g2, δ)
    @test integral_1 ≈ 0.0
    @test integral_2 ≈ 0.0
    integral_3 = _integrate_gap_numeric(ξ, g1, δ)
    @test isapprox(integral_1, integral_3, atol = 1e-7)
end

function test_integral_in_out_of_contact()
    ξ = SVector{2, Float64}(0.2, 0.9)
    g1 = SVector{2, Float64}(0.1, 0.7)
    g2 = SVector{2, Float64}(0.7, 0.1)
    δ = 0.5
    integral_1 = Cthonios.integrate_gap(ξ, g1, δ)
    integral_2 = Cthonios.integrate_gap(ξ, g2, δ)
    @test integral_1 ≈ integral_2
    integral_3 = Cthonios.integrate_gap(ξ, g1, δ)
    @test isapprox(integral_1, integral_3, atol = 1e-7)
end

@testset "ContactPair" begin
    test_contact_pair_constructor()
end

@testset "ContactSurface" begin
    test_contact_surface_constructor_and_show()
end

@testset "Nearest neighbor search" begin
    test_neareast_neighbor_search_single_quad()
    test_neareast_neighbor_facet_field_single_quad()
    for mn in [1, 2, 3, 4, 5]
        test_neareast_neighbor_search_two_block_tri3_mesh(mn)
    end
end

@testset "Integral penalty contact tests" begin
    test_integral_in_overlap()
    test_integral_both_in_contact()
    test_integral_both_equal_in_contact()
    test_integral_both_out_of_contact()
    test_integral_both_equal_out_of_contact()
    test_integral_in_out_of_contact()
end

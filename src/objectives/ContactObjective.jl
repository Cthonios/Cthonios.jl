struct PenaltyContactObjectiveCache{
    F, C,
    A, O, RT, RV, 
    P <: AbstractObjectiveCache{A, O, RT, RV}
} <: AbstractObjectiveCache{A, O, RT, RV}
    contact_formulation::F
    contact_pairs::C
    physics_cache::P
end

FiniteElementContainers.create_unknowns(o::PenaltyContactObjectiveCache) = 
create_unknowns(o.physics_cache)

assembler(o::PenaltyContactObjectiveCache) = assembler(o.physics_cache)

function gradient(o::PenaltyContactObjectiveCache, U, p)
    g = gradient(o.physics_cache, U, p)
    c = create_field(o.physics_cache)
    assemble_contact_vector!(c, o.contact_formulation, PenaltyContact(), o.contact_pairs, p.h1_coords, p.h1_field)
    return g .+ 1e3 * c.data
end

function hessian(o::PenaltyContactObjectiveCache, U, p; symmetric = true)
    # TODO eventually tack on contact contribution
    return hessian(o.physics_cache, U, p; symmetric = symmetric)
end 

function hvp(o::PenaltyContactObjectiveCache, U, V, p)
    return hvp(o.physics_cache, U, V, p)
end

function value(o::PenaltyContactObjectiveCache, U, p)
    e = value(o.physics_cache, U, p)
    c = zeros(1)
    assemble_contact_scalar!(c, o.contact_formulation, PenaltyContact(), o.contact_pairs, p.h1_coords, p.h1_field)
    return e + 1e3 * c[1]
end

function initialize!(o::PenaltyContactObjectiveCache, U, p)
    initialize!(o.physics_cache, U, p)
    return nothing
end

function step!(solver, o::PenaltyContactObjectiveCache, U, p; verbose = true)
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p)
    update_nearest_neighbors!(o.contact_pairs, p.h1_coords, U)
    _step_begin_banner(o.physics_cache, p; verbose = verbose)
    # do solve TODO
    solve!(solver, U.data, p)
    update_field_dirichlet_bcs!(U, p.dirichlet_bcs)

    # update values at end of step
    gradient(o, U, p)
    value(o, U, p)
    # update gradient/vale at end of step TODO

end

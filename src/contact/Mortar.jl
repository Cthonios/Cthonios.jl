Base.@kwdef struct MortarContact{RT} <: AbstractContactEnforcement
    max_overlap_dist::RT
    rel_smoothing_size::RT = 1.e-1
end

function assemble_contact_matrix!(
    m::MortarContact,
    ::PenaltyContact,
    cache::ContactPairCache,
    X, U
)
    side_a = cache.side_a
    side_b = cache.side_b
    interactions = cache.interactions

    energy = zero(eltype(X))
    for b in axes(side_b.side_nodes, 2)
        for a in axes(interactions, 1)
            index_a = interactions[a, b]
            invalid_index = index_a == -1

            if invalid_index
                # energy = zero(eltype(X))
                # do nothing for this case
                @warn "Have an invalid index"
            else
                X_a = facet_field(side_a, X, a)
                X_b = facet_field(side_b, X, b)
                U_a = facet_field(side_a, U, a)
                U_b = facet_field(side_b, U, b)

                display(facet_normal(X_a))
                display(facet_normal(X_b))
                # energy += contact_energy(X_a, X_b, U_a, U_b, m.max_overlap_dist, m.rel_smoothing_size)
                # K_d = contact_stiffness_driver(X_a, X_b, U_a, U_b, m.max_overlap_dist, m.rel_smoothing_size)
                K_f = contact_stiffness_follow(X_a, X_b, U_a, U_b, m.max_overlap_dist, m.rel_smoothing_size)
                # display(K_d)
                display(K_f)
            end
        end
    end
    return energy
end

function assemble_contact_scalar!(
    energy,
    m::MortarContact,
    type::PenaltyContact,
    cache::ContactPairCache,
    X, U
)
    block_energy = zero(eltype(X))
    side_a = cache.side_a
    side_b = cache.side_b
    interactions = cache.interactions
    for b in axes(side_b.side_nodes, 2)
        for a in axes(interactions, 1)
            index_a = interactions[a, b]
            side_a.side_nodes[:, index_a]
            invalid_index = index_a == -1

            if invalid_index
                # energy = zero(eltype(X))
                # do nothing for this case
                # @warn "Have an invalid index"
            else
                X_a = facet_field(side_a, X, index_a)
                X_b = facet_field(side_b, X, b)
                U_a = facet_field(side_a, U, index_a)
                U_b = facet_field(side_b, U, b)
                block_energy += contact_energy(type, X_a, X_b, U_a, U_b, m.max_overlap_dist, m.rel_smoothing_size)
            end
            # @assert false
        end
    end
    # fill!(contact_energy, energy)
    energy .+= block_energy
    return nothing
end

function assemble_contact_vector!(
    vector,
    m::MortarContact,
    type::PenaltyContact,
    cache::ContactPairCache,
    X, U
)
    side_a = cache.side_a
    side_b = cache.side_b
    interactions = cache.interactions
    n_dofs = size(U, 1)
    energy = zero(eltype(X))
    for b in axes(side_b.side_nodes, 2)
        for a in axes(interactions, 1)
            index_a = interactions[a, b]
            invalid_index = index_a == -1

            if invalid_index
                # energy = zero(eltype(X))
                # do nothing for this case
                @warn "Have an invalid index"
            else
                X_a = facet_field(side_a, X, index_a)
                X_b = facet_field(side_b, X, b)
                U_a = facet_field(side_a, U, index_a)
                U_b = facet_field(side_b, U, b)
                R_d = contact_residual_side_a(type, X_a, X_b, U_a, U_b, m.max_overlap_dist, m.rel_smoothing_size)
                R_f = contact_residual_side_b(type, X_a, X_b, U_a, U_b, m.max_overlap_dist, m.rel_smoothing_size)

                # now assemble the two contributions
                local_dof = 1
                for node in axes(side_a.side_nodes, 1)
                    for dof in 1:n_dofs
                        vector[dof, side_a.side_nodes[node, index_a]] = R_d.data[local_dof]
                        local_dof = local_dof + 1
                    end
                end

                local_dof = 1
                for node in axes(side_b.side_nodes, 1)
                    for dof in 1:n_dofs
                        vector[dof, side_b.side_nodes[node, b]] = R_f.data[local_dof]
                        local_dof = local_dof + 1
                    end
                end
            end
        end
    end
    return energy
end

function contact_residual_side_a(
    type::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.gradient(
        z -> contact_energy(
            type,
            X_a, X_b, z, U_b, 
            max_overlap_dist, rel_smoothing_size
        ), U_a
    )
end

function contact_residual_side_b(
    type::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.gradient(
        z -> contact_energy(
            type,
            X_a, X_b, U_a, z, 
            max_overlap_dist, rel_smoothing_size
        ), U_b
    )
end

function contact_stiffness_side_a_side_a(
    type::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.hessian(
        z -> contact_energy(
            type,
            X_a, X_b, z, U_b, 
            max_overlap_dist, rel_smoothing_size
        ), U_a
    )
end

function contact_stiffness_side_b_side_b(
    type::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.hessian(
        z -> contact_energy(
            type,
            X_a, X_b, U_a, z, 
            max_overlap_dist, rel_smoothing_size
        ), U_b
    )
end

function contact_stiffness_side_a_side_b(
    type::AbstractContactFormulation,
    X_a, X_b, U_a, U_b,
    max_overlap_dist, 
    rel_smoothing_size
)
    return ForwardDiff.jacobian(
        z ->  contact_residual_side_a(
                type, 
                X_a, X_b, U_a, z, 
                max_overlap_dist, rel_smoothing_size
            ), U_b
        )
end

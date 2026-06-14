mutable struct EigenObjective{RT, RV, RM <: AbstractMatrix{RT}, A} <: AbstractSolidMechanicsObjective{RT, RV, A}
    assembler::A
    eigen_vals::RV
    eigen_vecs::RM
    timer::TimerOutput
end

function EigenObjective(assembler, timer::TimerOutput = TimerOutput())
    eigen_vals = zeros(0)
    eigen_vecs = zeros(length(assembler.residual_storage), 0)
    timer = TimerOutput()
    return EigenObjective(assembler, eigen_vals, eigen_vecs, timer)
end

# TODO we can probably just remove times all together
function initialize!(::EigenObjective, u, p)
    # p.times.time_current = zero(typeof(p.times.time_current))
    return nothing
end

function postprocess!(
    pp, output_settings, n, objective::EigenObjective, u, p
)
    nat_circ_freqs = sqrt.(1 ./ objective.eigen_vals)
    for n in axes(nat_circ_freqs, 1)
        eig_vec_unknowns = objective.eigen_vecs[:, n]
        u_field = create_field(objective)
        update_field_unknowns!(u_field, objective.assembler.dof, u)
        update_field_dirichlet_bcs!(u_field, p.dirichlet_bcs)
        unknown_dofs = objective.assembler.dof.unknown_dofs
        FEC.fec_foreach(eig_vec_unknowns) do dof
            u_field[unknown_dofs[dof]] += eig_vec_unknowns[dof]
        end

        write_times(pp, n, nat_circ_freqs[n])
        if size(p.field, 1) == 2
            write_field(pp, n, ("displ_x", "displ_y"), u_field)
        elseif size(p.field, 1) == 3
            write_field(pp, n, ("displ_x", "displ_y", "displ_z"), u_field)
        end
    end
end

function run!(
    solver, objective::EigenObjective, u, p,
    mesh, output_file,
    output_settings = default_output_settings(objective),
    verbose::Bool = true
)
    pp = initialize_postprocessor(mesh, output_file, objective, output_settings)
    initialize!(objective, u, p)
    step!(solver, objective, u, p; verbose = verbose)
    try
        postprocess!(pp, output_settings, 1, objective, u, p)
    finally
        close(pp)
    end
    return nothing
end

function step!(solver, o::EigenObjective, u, p; verbose = false)
    solve!(solver, o, u, p)
    return nothing
end

function hessian(o::EigenObjective, u, p)
    assemble_stiffness!(assembler(o), stiffness!, u, p)
    return stiffness(assembler(o))
end

function mass_matrix(o::EigenObjective, u, p)
    assemble_mass!(assembler(o), mass!, u, p)
    return FEC.mass(assembler(o))
end

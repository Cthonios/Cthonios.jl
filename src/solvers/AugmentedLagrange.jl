Base.@kwdef struct AugmentedLagrangeSolverSettings
    penalty_scaling::Float64 = 4.
    target_constraint_decrease_factor::Float64 = 0.75
    relative_gmres_tol::Float64 = 2.e-2
    max_gmres_iters::Int = 100
    use_second_order_update::Bool = true
    use_newton_only::Bool = false
    num_initial_low_order_iterations::Int = 3
    inverse_ncp_hessian_bound::Float64 = 1.e-2
    max_al_iters::Int = 100
    tol::Float64 = 1.e-8
end

struct AugmentedLagrangeSolver
    settings::AugmentedLagrangeSolverSettings
end

function AugmentedLagrangeSolver(objective_cache, p; kwargs...)
    settings = AugmentedLagrangeSolverSettings(; kwargs...)
    return AugmentedLagrangeSolver(settings)
end

function linear_update()

end

function solve!(solver::AugmentedLagrangeSolver, Uu, p, λ, κ)
    # TODO add warm start

    # TODO add update preconditioner

    huge_val = 1.e64
    ncp_error = huge_val * ones(length(λ))
    error_norm = huge_val

    initial_tol_scaling = 100.0
    tol_ramp_down_iters = 3

    for it in 1:solver.settings.max_al_iters
        # TODO finish this up... we'll need a way
        # to dynamically change equation solver settings
        # if it < tol_ramp_down_iters
        #     tol_scaling = initial_tol_scaling^((tol_ramp_down_iters - it) / tol_ramp_down_iters)

        # end

        update_precond = false
        if (solver.settings.use_second_order_update && it >= solver.settings.num_initial_low_order_iterations) ||
            solver.settings.use_newton_only

            # do linear update
            # TODO finish this
            dx, dl, solve_success = linear_update(objective_cache, Uu, p, λ, κ)

            if !solve_success
                update_precond = true
            end

            λ_save = copy(λ)

            # TODO move to settings
            max_line_search_iters = 10
            for linesearch in 1:max_line_search_iters
                if linesearch == max_line_search_iters
                    update_precond = true
                end

                y = Uu + dx
                λ .+= dl

                # TODO need total residual call here
                trial_error_norm = norm()

                if norm(trial_error_norm < error_norm)
                    @info "Total error after 2nd order update = $trial_error_norm"
                    error_norm = trial_error_norm
                    copyto!(Uu, y)
                else
                    @info "Total error after 2nd order update = $trial_error_norm, no improvement"
                    copyto!(λ, λ_save)
                    dx .*= 0.2
                    dl .*= 0.2
                end
            end
        end

        if update_precond
            # TODO make call to update precond
        end

        if !solver.settings.use_newton_only

        end
    end
end

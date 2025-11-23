Base.@kwdef struct NewtonSolverSettings
    max_iters::Int = 50
    rel_tol::Float64 = 1e-5
    abs_tol::Float64 = 1e-6
    verbose::Bool = true
end

struct NewtonSolver{O}
    objective_cache::O
    settings::NewtonSolverSettings
end

function NewtonSolver(objective_cache; kwargs...)
    return NewtonSolver(objective_cache, NewtonSolverSettings(; kwargs...))
end 

function solve!(solver::NewtonSolver, Uu, p)
    objective_cache = solver.objective_cache
    # R = gradient(objective_cache, Uu, p)
    # res_norm0 = norm(R)

    # if res_norm0 == 0.0
    #     res_norm0 = 1.
    # end
    res_norm0 = Inf
    n = 1

    if solver.settings.verbose
        @info "=============================================================================================================================="
        @info "Iteration |R|          |R| / |R0|   |dUu|        Reference"
        @info "=============================================================================================================================="
    end

    while n <= solver.settings.max_iters
        R = gradient(objective_cache, Uu, p)

        if n == 1
            res_norm0 = norm(R)
        end

        K = hessian(objective_cache, Uu, p)
        dUu = -K \ R
        # @show norm(R) norm(R) / res_norm0 norm(dUu) 
        # copyto!(solver.solution, -K \ R)
        res_norm = norm(R)
        rel_res_norm = res_norm / res_norm0

        if solver.settings.verbose
            str = @sprintf "%8d  %1.6e %1.6e %1.6e %1.6e" n res_norm rel_res_norm norm(dUu) res_norm0
            @info str
        end

        if rel_res_norm < solver.settings.rel_tol ||
            res_norm < 1.e-8
            Uu .+= dUu
            if solver.settings.verbose
                @info "Converged"
            end
            return nothing
        end
        # objective.solution .+= dUu
        Uu .+= dUu

        # if norm(R) 
        n = n + 1
    end

    error("Newton solver failed to converge")
end

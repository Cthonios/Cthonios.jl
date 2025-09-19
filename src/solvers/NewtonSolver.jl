Base.@kwdef struct NewtonSolverSettings
    max_iters::Int = 50
    rel_tol::Float64 = 1e-5
    abs_tol::Float64 = 1e-6
end

struct NewtonSolver{O, S}
    objective_cache::O
    settings::S
end

function NewtonSolver(objective_cache)
    return NewtonSolver(objective_cache, NewtonSolverSettings())
end 

function solve!(solver::NewtonSolver, Uu, p)
    objective_cache = solver.objective_cache
    R = gradient(objective_cache, Uu, p)
    res_norm0 = norm(R)

    if res_norm0 == 0.0
        res_norm0 = 1.
    end
    n = 1

    @info "=============================================================================================================================="
    @info "Iteration |R|          |R| / |R0|   |dUu|        Reference"
    @info "=============================================================================================================================="

    while n <= solver.settings.max_iters
        R = gradient(objective_cache, Uu, p)
        K = hessian(objective_cache, Uu, p)
        dUu = -K \ R
        # @show norm(R) norm(R) / res_norm0 norm(dUu) 
        # copyto!(solver.solution, -K \ R)
        res_norm = norm(R)
        rel_res_norm = res_norm / res_norm0
        str = @sprintf "%8d  %1.6e %1.6e %1.6e %1.6e" n res_norm rel_res_norm norm(dUu) res_norm0
        @info str
        if rel_res_norm < solver.settings.rel_tol ||
            res_norm < 1.e-8
            Uu .+= dUu
            @show "Converged"
            return nothing
        end
        # objective.solution .+= dUu
        Uu .+= dUu

        # if norm(R) 
        n = n + 1
    end

    error("Newton solver failed to converge")
    return nothing
end

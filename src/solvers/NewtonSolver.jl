Base.@kwdef struct NewtonSolverSettings
    max_iters::Int = 10
    rel_tol::Float64 = 1e-8
    abs_tol::Float64 = 1e-8
end

struct NewtonSolver{O, S}
    objective::O
    settings::S
end

function NewtonSolver(objective)
    return NewtonSolver(objective, NewtonSolverSettings())
end 

function solve!(solver::NewtonSolver, Uu, p)
    objective = solver.objective
    R = gradient(objective, Uu, p)
    res_norm0 = norm(R)

    if res_norm0 == 0.0
        res_norm0 = 1.
    end
    n = 1

    @info "=============================================================================================================================="
    @info "Iteration |R|          |R| / |R0|   |dUu|        Reference"
    @info "=============================================================================================================================="

    while n <= solver.settings.max_iters
        R = gradient(objective, Uu, p)
        K = hessian(objective, Uu, p)
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
            break
        end
        # objective.solution .+= dUu
        Uu .+= dUu

        # if norm(R) 
        n = n + 1
    end
    return nothing
end

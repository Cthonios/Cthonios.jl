Base.@kwdef struct NewtonSolverSettings
    max_iters::Int = 50
    rel_tol::Float64 = 1e-8
    abs_tol::Float64 = 1e-8
    verbose::Bool = true
end

struct NewtonSolver{L}
    linear_solver::L
    settings::NewtonSolverSettings
end

function NewtonSolver(objective, u, p; kwargs...)
    linear_solver = KrylovSolver(Val{:cg}(), objective, u, p)
    return NewtonSolver(linear_solver, NewtonSolverSettings(; kwargs...))
end 

function solve!(solver::NewtonSolver, objective, u, p, P)
    res_norm0 = one(eltype(u))
    n = 1

    if solver.settings.verbose
        @info "=============================================================================================================================="
        @info "Iteration |R|          |R| / |R0|   |dUu|        Reference"
        @info "=============================================================================================================================="
    end

    while n <= solver.settings.max_iters
        R = gradient(objective, u, p)

        if n == 1
            norm_R = norm(R)
            if norm_R != 0.
                res_norm0 = norm(R)
            end
        end

        # K = hessian(objective_cache, u, p)
        # dUu = -K \ R
        Δu = Cthonios.solve!(solver.linear_solver, objective, u, p, P)
        norm_du = norm(Δu)
        # @show norm(R) norm(R) / res_norm0 norm(dUu) 
        # copyto!(solver.solution, -K \ R)
        res_norm = norm(R)
        rel_res_norm = res_norm / res_norm0

        if solver.settings.verbose
            str = @sprintf "%8d  %1.6e %1.6e %1.6e %1.6e" n res_norm rel_res_norm norm_du res_norm0
            @info str
        end

        if rel_res_norm < solver.settings.rel_tol ||
            res_norm < 1.e-8
            # Uu .+= dUu
            if solver.settings.verbose
                @info "Converged"
            end
            return nothing
        end
        # objective.solution .+= dUu
        # Uu .+= dUu
        u .-= Δu

        # if norm(R) 
        n = n + 1
    end

    error("Newton solver failed to converge")
end

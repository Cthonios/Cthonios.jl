Base.@kwdef struct ExplicitSolverSettings
    verbose::Bool = true
end

struct ExplicitSolver{O}
    objective_cache::O
    settings::ExplicitSolverSettings
end

function ExplicitSolver(objective_cache, p)
    return ExplicitSolver(objective_cache, ExplicitSolverSettings())
end

function solve!(solver::ExplicitSolver, u, p)
    o = solver.objective_cache
    @timeit o.timer "gradient" begin
        R_eff = gradient(o, u, p)
    end
    m = lumped_mass(o, u, p)
    @timeit o.timer "solve" begin
        @. o.a = -(R_eff / m)
    end
    return nothing
end

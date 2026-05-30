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

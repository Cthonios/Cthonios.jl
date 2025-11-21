Base.@kwdef struct EigenSolverSettings
    n_eigs::Int
    max_iters::Int = 1000
end

struct EigenSolver{O, S}
    objective_cache::O
    settings::S
end

function EigenSolver(objective_cache, n_eigs::Int)
    return EigenSolver(objective_cache, EigenSolverSettings(n_eigs=n_eigs))
end

function solve!(solver::EigenSolver)
    objective_cache = solver.objective_cache
    Uu, p = objective_cache.solution, objective_cache.parameters
    K = hessian(objective_cache, Uu, p).data
    M = mass_matrix(objective_cache, Uu, p).data
    vals, vecs = eigs(M, K, nev=solver.settings.n_eigs, which=:LM)
end

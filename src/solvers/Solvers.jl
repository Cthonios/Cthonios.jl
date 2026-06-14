include("LinearSolvers.jl")
include("Preconditioners.jl")
include("Predictors.jl")

include("NewtonSolver.jl")
include("TrustRegionSolver.jl")

include("AugmentedLagrange.jl")

Base.@kwdef struct EigenSolverSettings
    n_eigs::Int
    max_iters::Int = 1000
end

struct EigenSolver
    settings::EigenSolverSettings
    timer::TimerOutput
end

function EigenSolver(objective, u, p, n_eigs; timer::TimerOutput = TimerOutput())
    return EigenSolver(
        EigenSolverSettings(n_eigs = n_eigs), timer
    )
end

function solve!(solver::EigenSolver, objective, u, p)
    n_eigs = solver.settings.n_eigs
    K = hessian(objective, u, p)
    M = mass_matrix(objective, u, p)
    vals, vecs = eigs(M, K, nev = n_eigs, which = :LM)
    objective.eigen_vals = vals
    objective.eigen_vecs = vecs
    return nothing
end

Base.@kwdef struct ExplicitSolverSettings
    verbose::Bool = true
end

struct ExplicitSolver
    settings::ExplicitSolverSettings
end

function ExplicitSolver(objective, u, p)
    return ExplicitSolver(ExplicitSolverSettings())
end

Base.@kwdef struct ImplicitSolverSettings
    verbose::Bool = true
end

struct ImplicitSolver{
    N,
    P1 <: AbstractPreconditioner,
    P2 <: AbstractPredictor,
}
    nonlinear_solver::N
    preconditioner::P1
    predictor::P2
    settings::ImplicitSolverSettings
    timer::TimerOutput

    function ImplicitSolver(
        nonlinear_solver, preconditioner, predictor;
        settings = ImplicitSolverSettings(),
        timer = TimerOutput()
    )
        new{typeof(nonlinear_solver), typeof(preconditioner), typeof(predictor)}(
            nonlinear_solver, preconditioner, predictor, settings, timer
        )
    end
end


function solve!(solver::ImplicitSolver, objective, u, p)
    P = solver.preconditioner
    @info "ImplicitSolver"
    @info "  Predictor"
    solve!(solver.predictor, objective, u, p, P)
    @info "  Nonlinear solve"
    solve!(solver.nonlinear_solver, objective, u, p, P)
    return nothing
end

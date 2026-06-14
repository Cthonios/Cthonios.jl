struct KrylovSolverSettings <: FEC.AbstractLinearSolverSettings
end

struct KrylovSolver{
    # P,
    W <: Krylov.KrylovWorkspace
}
    # preconditioner::P
    settings::KrylovSolverSettings
    timer::TimerOutput
    workspace::W

    function KrylovSolver(
        solver_type::Val{Sym},
        objective::AbstractObjective{RT, RV, A}, u, p;
        # preconditioner;
        settings::KrylovSolverSettings = KrylovSolverSettings(),
        timer::TimerOutput = TimerOutput()
    ) where {Sym, RT, RV, A}
        g = gradient(objective, u, p)
        H = hessian(objective, u, p)
        workspace = krylov_workspace(solver_type, H, g)

        # TODO probably want a better way to setup preconditioner
        # new{typeof(preconditioner), typeof(workspace)}(
        new{typeof(workspace)}(
            # preconditioner, settings, timer, workspace
            settings, timer, workspace
        )
    end

    function KrylovSolver(
        solver_type::Val{Sym},
        H::AbstractMatrix,
        g::AbstractVector;
        # preconditioner;
        settings::KrylovSolverSettings = KrylovSolverSettings(),
        timer::TimerOutput = TimerOutput()
    ) where Sym
        workspace = Krylov.krylov_workspace(solver_type, H, g)
        # new{typeof(preconditioner), typeof(workspace)}(
        #     preconditioner, settings, timer, workspace
        new{typeof(workspace)}(
            settings, timer, workspace
        )
    end
end

function solve!(solver::KrylovSolver, H::AbstractMatrix, g, P)
    @timeit solver.timer "KrylovSolver - solve!" begin
        @timeit solver.timer "krylov_solve!" begin
            Krylov.krylov_solve!(solver.workspace, H, g, M = P, ldiv = true)
        end
        Δu, stats = Krylov.results(solver.workspace)
        return Δu
    end
end

function solve!(solver::KrylovSolver, objective::AbstractObjective, u, p, P)
    @timeit solver.timer "KrylovSolver - solve!" begin
        g = gradient(objective, u, p)
        H = hessian(objective, u, p)
        @timeit solver.timer "krylov_solve!" begin
            Krylov.krylov_solve!(solver.workspace, H, g, M = P, ldiv = true)
        end
        # TODO we can add some logging
        Δu, stats = Krylov.results(solver.workspace)
    return Δu
    end
end

abstract type AbstractSolver end

"""
$(TYPEDSIGNATURES)
Creates a set of unknowns for the nonlinear solver
"""
function create_unknowns(solver::AbstractSolver)
  return create_unknowns(solver.objective.domain)
end

include("eigen_solvers/EigenSolvers.jl")
include("linear_solvers/LinearSolvers.jl")
include("nonlinear_solvers/NonlinearSolvers.jl")
# include("WarmStart.jl")

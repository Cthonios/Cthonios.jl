include("linear_solvers/LinearSolvers.jl")

abstract type NonlinearSolverSettings end
abstract type NonlinearSolver{S, L} end

function logger end # To be defined for each solver
function solve! end # TO be defined for each solver
function update_unknown_dofs! end # To be defined for each solver

function update_unknown_dofs!(solver::NonlinearSolver, d::QuasiStaticDomain)
  # update the dofs using the linear solver method
  update_unknown_dofs!(solver.linear_solver, d)
end

include("NewtonSolver.jl")
include("TrustRegionSolver.jl")

function setup_nonlinear_solver(
  inputs::D, domain::QuasiStaticDomain, backend
) where {D <: Dict{Symbol, Any}}
  type = eval(Meta.parse(inputs[:type]))
  return type(inputs, domain, backend)
end

include("LinearSolvers.jl")

abstract type NonlinearSolverSettings end
abstract type NonlinearSolver{S, L, V} end

# for changing size after bc changes
function Base.resize!(solver::NonlinearSolver, domain::QuasiStaticDomain)
  resize!(solver.Uu, length(domain.dof.unknown_dofs))
  resize!(solver.ΔUu, length(domain.dof.unknown_dofs))
end
function update_increment!(
  solver::NonlinearSolver, ΔUu::V
) where V <: AbstractVector
  solver.ΔUu .= ΔUu
end
function logger end # TO be defined for each solver
function solve! end # TO be defined for each solver

include("NewtonSolver.jl")
include("TrustRegionSolver.jl")

function setup_nonlinear_solver(
  inputs::D, domain::QuasiStaticDomain
) where {D <: Dict{Symbol, Any}}
  type = eval(Meta.parse(inputs[:type]))
  return type(inputs, domain)
end

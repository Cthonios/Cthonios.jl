include("linear_solvers/LinearSolvers.jl")

abstract type NonlinearSolverSettings end
abstract type NonlinearSolver{S, L, V} end

# for changing size after bc changes
function Base.resize!(solver::NonlinearSolver, domain::QuasiStaticDomain)
  resize!(solver.Uu, length(domain.dof.unknown_dofs))
  resize!(solver.ΔUu, length(domain.dof.unknown_dofs))
end
# function update_increment!(
#   solver::NonlinearSolver, ΔUu::V
# ) where V <: AbstractVector
#   solver.ΔUu .= ΔUu
# end
function logger end # TO be defined for each solver
function solve! end # TO be defined for each solver

function update_unknown_dofs!(solver::NonlinearSolver, d::QuasiStaticDomain)
  # update the dofs
  FiniteElementContainers.update_unknown_dofs!(d.dof, d.bc_dofs)
  FiniteElementContainers.update_unknown_dofs!(solver.linear_solver.assembler, d.dof, map(x -> x.fspace, d.sections), d.bc_dofs)

  # now resize the caches
  resize!(solver.Uu, length(d.dof.unknown_dofs))
  resize!(solver.ΔUu, length(d.dof.unknown_dofs))
end

include("NewtonSolver.jl")
include("TrustRegionSolver.jl")

function setup_nonlinear_solver(
  inputs::D, domain::QuasiStaticDomain, backend
) where {D <: Dict{Symbol, Any}}
  type = eval(Meta.parse(inputs[:type]))
  return type(inputs, domain, backend)
end

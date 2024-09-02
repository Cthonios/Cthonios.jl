"""
$(TYPEDEF)
a nonlinear solver needs to define the following
methods

1. check_convergence - returns a bool 
2. logger - @info information
3. step! - step for the solver

it needs the following types
1. a linear solver
2. an objective
3. an unknown vector
4. an int called max_iter
"""
abstract type AbstractNonlinearSolver{L, O, U} end

"""
$(TYPEDSIGNATURES)
Creates a set of unknowns for the nonlinear solver
"""
function create_unknowns(solver::AbstractNonlinearSolver)
  return create_unknowns(solver.objective.domain)
end

"""
$(TYPEDSIGNATURES)
Computes the residual and in-place updates
the appropriate storage in the linear solver.
This method assumes solver.U has already been updated
so ``update_fields!`` should be called prior to calling
this.
"""
function residual!(solver::AbstractNonlinearSolver, p)
  residual!(solver.linear_solver, solver.objective, p)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function solve!(solver::AbstractNonlinearSolver, Uu, p)
  # calculate initial residual for logging and 
  # convergence check purposes
  update_fields!(solver, Uu)
  residual!(solver, p)
  norm_R0 = residual_norm(solver.linear_solver)

  # loop over nonlinear iterations
  for n in 1:solver.max_iter
    step!(solver, Uu, p)
    Uu .= Uu .- solver.ΔUu
    logger(solver, n, norm_R0)

    if check_convergence(solver, norm_R0)
      return nothing
    end
  end
  @info "Maximum iterations reached"
  return nothing
end

"""
$(TYPEDSIGNATURES)
Computes the tangents and in-place updates
the appropriate storage in the linear solver.
This method assumes solver.U has already been updated
so ``update_fields!`` should be called prior to calling
this.

This method may not be defined for all linear solvers.
"""
function tangent!(solver::AbstractNonlinearSolver, p)
  tangent!(solver.linear_solver, solver.objective, p)
  return nothing
end

"""
$(TYPEDSIGNATURES)
Wraps update_fields! for the linear solver to
update the unknowns in solver.linear_solver.U
with the values in Uu
"""
function update_fields!(solver::AbstractNonlinearSolver, Uu)
  update_fields!(solver.linear_solver, solver.objective, Uu)
  return nothing
end

# nonlinear solvers
include("NewtonSolver.jl")

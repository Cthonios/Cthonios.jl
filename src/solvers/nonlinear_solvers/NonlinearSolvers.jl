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
abstract type AbstractNonlinearSolver{L, O, U, T} end
timer(s::AbstractNonlinearSolver) = s.timer

"""
$(TYPEDSIGNATURES)
Creates a set of unknowns for the nonlinear solver
"""
function create_unknowns(solver::AbstractNonlinearSolver)
  return create_unknowns(solver.objective.domain)
end

function objective(solver::AbstractNonlinearSolver, Uu, p)
  return objective!(solver.o, solver.objective, Uu, p)
end

function gradient(solver::AbstractNonlinearSolver, Uu, p)
  return gradient!(solver.g, solver.objective, Uu, p)
end

function hvp(solver::AbstractNonlinearSolver, Uu, p, Vv)
  return hvp!(solver.Hv, solver.objective, Uu, p, Vv)
end

"""
$(TYPEDSIGNATURES)
Generic method to fall back on if step! is defined
"""
function solve!(solver::AbstractNonlinearSolver, Uu, p)
  @timeit timer(solver) "AbstractNonlinearSolver - solve!" begin
    gradient!(solver.linear_solver, solver.objective, Uu, p)
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
  end
  return nothing
end

# nonlinear solvers
include("NewtonSolver.jl")
include("TrustRegionSolver.jl")

# exports
export NewtonSolver
export TrustRegionSolver

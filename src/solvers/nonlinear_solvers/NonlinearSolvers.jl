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
abstract type AbstractNonlinearSolver{L, O, U, W, T} <: AbstractSolver end
timer(s::T) where T <: AbstractNonlinearSolver = s.timer

# """
# $(TYPEDSIGNATURES)
# Generic method to fall back on if step! is defined
# """
# function solve!(solver::AbstractNonlinearSolver, Uu, p)
#   @timeit timer(solver) "AbstractNonlinearSolver - solve!" begin
#     gradient!(solver.linear_solver, solver.objective, Uu, p)
#     norm_R0 = residual_norm(solver.linear_solver, solver.objective, Uu, p)

#     # loop over nonlinear iterations
#     for n in 1:solver.max_iter
#       step!(solver, Uu, p)
#       Uu .= Uu .- solver.ΔUu
#       logger(solver, Uu, p, n, norm_R0)

#       if check_convergence(solver, Uu, p, norm_R0)
#         return nothing
#       end
#     end
#     @info "Maximum iterations reached"
#   end
#   return nothing
# end

# nonlinear solvers
include("WarmStart.jl")

# include("NewtonSolver.jl")
include("TrustRegionSolver.jl")

# exports
# export NewtonSolver
export TrustRegionSolver

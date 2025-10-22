"""
$(TYPEDEF)
A linear solver simply needs to define the
following methods
1. residual_norm
2. solve!
solve! method which has the following signature
solve!(ΔUu, solver::DirectSolver, obj, Uu, p)
where ΔUu is the increment for a nonlinear solver
solver is the solver, obj is the objective,
Uu is the current guess of the solution, and
p is the set of parameters
"""
abstract type AbstractLinearSolver end

# linear solvers
# include("DirectSolver.jl")
# include("KrylovSolver.jl")
include("Preconditioners.jl")

# preconditioners
export CholeskyPreconditioner
export LimitedLDLPreconditioner
export LDLPreconditioner

# solvers
# export DirectSolver
# export KrylovSolver

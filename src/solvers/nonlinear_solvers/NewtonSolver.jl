"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct NewtonSolver{L, O, U, T} <: AbstractNonlinearSolver{L, O, U, T}
  linear_solver::L
  objective::O
  ΔUu::U
  timer::T
  max_iter::Int
  abs_tol::Float64
  rel_tol::Float64
  use_warm_start::Bool
end

"""
$(TYPEDSIGNATURES)
"""
function NewtonSolver(objective::Objective, p, linear_solver_type, timer)
  @timeit timer "NewtonSolver - setup" begin
    linear_solver = linear_solver_type(objective, p, timer)
    ΔUu = create_unknowns(objective.domain)
  end
  return NewtonSolver(
    linear_solver, objective, ΔUu, timer,
    100, 1e-8, 1e-10, false
  )
end

"""
$(TYPEDSIGNATURES)
"""
function NewtonSolver(inputs::Dict{Symbol, Any}, objective::Objective, p, timer)
  @timeit timer "NewtonSolver - setup" begin
    linear_solver_inputs = inputs[Symbol("linear solver")]
    linear_solver = eval(Symbol(linear_solver_inputs[:type]))(linear_solver_inputs[:parameters], objective, p, timer)
    ΔUu = create_unknowns(objective.domain)
  end
  return NewtonSolver(
    linear_solver, objective, ΔUu, timer,
    100, 1e-8, 1e-10, false
  )
end

"""
$(TYPEDSIGNATURES)
"""
function check_convergence(solver::NewtonSolver, Uu, p, R0_norm)
  R_norm = residual_norm(solver.linear_solver, solver.objective, Uu, p)
  U_norm = norm(solver.ΔUu)

  if R_norm / R0_norm < solver.rel_tol || 
     U_norm < solver.rel_tol
    @info "Converged on relative tolerance"
    return true
  end

  if R_norm / R0_norm < solver.abs_tol || 
     U_norm < solver.abs_tol
    @info "Converged on absolute tolerance"
    return true
  end

  return false
end

"""
$(TYPEDSIGNATURES)
"""
function logger(solver::NewtonSolver, Uu, p, n::Int, norm_R0)
  norm_R = residual_norm(solver.linear_solver, solver.objective, Uu, p)
  norm_U = norm(solver.ΔUu)
  @info @sprintf "  Iteration %5i: ||R|| = %1.6e    ||R/R0|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R (norm_R / norm_R0) norm_U
end

"""
$(TYPEDSIGNATURES)
"""
function step!(solver::NewtonSolver, Uu, p)
  @timeit timer(solver) "NewtonSolver - step!" begin
    solve!(solver.ΔUu, solver.linear_solver, solver.objective, Uu, p)
  end
  return nothing
end

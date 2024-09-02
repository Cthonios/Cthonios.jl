"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct NewtonSolver{L, O, U} <: AbstractNonlinearSolver{L, O, U}
  linear_solver::L
  objective::O
  ΔUu::U
  max_iter::Int
  abs_tol::Float64
  rel_tol::Float64
end

"""
$(TYPEDSIGNATURES)
"""
function NewtonSolver(objective::Objective, linear_solver_type)
  linear_solver = linear_solver_type(objective.domain)
  ΔUu = create_unknowns(objective.domain)
  return NewtonSolver(
    linear_solver, objective, ΔUu,
    10, 1e-8, 1e-10
  )
end

"""
$(TYPEDSIGNATURES)
"""
function check_convergence(solver::NewtonSolver, R0_norm)
  R_norm = residual_norm(solver.linear_solver)
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
function logger(solver::NewtonSolver, n::Int, norm_R0)
  norm_R = residual_norm(solver.linear_solver)
  norm_U = norm(solver.ΔUu)
  @info @sprintf "  Iteration %5i: ||R|| = %1.6e    ||R/R0|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R (norm_R / norm_R0) norm_U
end

"""
$(TYPEDSIGNATURES)
"""
function step!(solver::NewtonSolver, Uu, p)
  update_fields!(solver, Uu)
  solve!(solver.ΔUu, solver.linear_solver, solver.objective, Uu, p)
  return nothing
end

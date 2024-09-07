"""
$(TYPEDFIELDS)
Direct linear solver
"""
struct DirectSolver{A, L} <: AbstractLinearSolver
  assembler::A
  linsolve::L
end

"""
$(TYPEDSIGNATURES)
Direct linear solver constructor from a Domain
"""
function DirectSolver(domain::Domain)
  asm = StaticAssembler(domain)
  update_unknown_dofs!(domain, asm)
  R = asm.residuals[domain.dof.unknown_dofs]
  K = SparseArrays.sparse!(asm) |> Symmetric
  assumptions = OperatorAssumptions(false)
  prob = LinearProblem(K, R; assumptions=assumptions)
  linsolve = init(prob)
  return DirectSolver(asm, linsolve)
end

function Base.show(io::IO, sol::DirectSolver)
  println(io, "DirectSolver:")
  println(io, "  Algorithm            = $(sol.linsolve.alg)")
  println(io, "  Assumptions          = $(sol.linsolve.assumptions)")
  println(io, "  Left preconditioner  = $(sol.linsolve.Pl)")
  println(io, "  Right preconditioner = $(sol.linsolve.Pr)")
end

"""
$(TYPEDSIGNATURES)
"""
function gradient!(solver::DirectSolver, obj::Objective, Uu, p::ObjectiveParameters)
  R = solver.assembler.residuals
  # annoying below
  update_fields!(p.U, obj.domain, Uu)
  R .= zero(eltype(R))
  domain_iterator!(R, obj.gradient, obj.domain, p.U, p.X)
  # set filled residual in linear solve init
  solver.linsolve.b = R[obj.domain.dof.unknown_dofs]
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function hessian!(solver::DirectSolver, obj, Uu, p)
  K = hessian!(solver.assembler, obj, Uu, p)
  solver.linsolve.A = K
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function residual_norm(solver::DirectSolver)
  return norm(solver.linsolve.b)
end

"""
$(TYPEDSIGNATURES)
"""
function solve!(ΔUu, solver::DirectSolver, obj, Uu, p)
  gradient!(solver, obj, Uu, p)
  hessian!(solver, obj, Uu, p)
  sol = LinearSolve.solve!(solver.linsolve)  
  ΔUu .= sol.u
  return nothing
end

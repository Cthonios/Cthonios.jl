"""
$(TYPEDFIELDS)
Direct linear solver
"""
struct DirectSolver{A, L, U1} <: AbstractLinearSolver
  assembler::A
  linsolve::L
  U::U1
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
  assumptions = OperatorAssumptions(true)
  prob = LinearProblem(K, R; assumptions=assumptions)
  linsolve = init(prob)
  U = create_fields(domain)
  return DirectSolver(asm, linsolve, U)
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
function residual!(solver::DirectSolver, obj, p)
  domain = obj.domain
  R = solver.assembler.residuals
  R .= zero(eltype(R))
  domain_iterator!(R, gradient, domain, solver.U, p)
  R_free = R[domain.dof.unknown_dofs]
  solver.linsolve.b = R_free
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
  update_fields!(solver.U, obj.domain, Uu)
  residual!(solver, obj, p)
  tangent!(solver, obj, p)
  sol = LinearSolve.solve!(solver.linsolve)  
  ΔUu .= sol.u
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function tangent!(solver::DirectSolver, obj, p)
  domain = obj.domain
  K = solver.assembler
  K.stiffnesses .= zero(eltype(K.stiffnesses))
  domain_iterator!(K, hessian, domain, solver.U, p)
  K = SparseArrays.sparse!(K) |> Symmetric
  solver.linsolve.A = K
  return nothing
end

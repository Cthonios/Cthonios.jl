"""
$(TYPEDFIELDS)
Direct linear solver
"""
struct DirectSolver{A, L, T} <: AbstractLinearSolver
  assembler::A
  linsolve::L
  timer::T
end
timer(s::DirectSolver) = s.timer

"""
$(TYPEDSIGNATURES)
Direct linear solver constructor from a Domain
"""
function DirectSolver(obj::Objective, p, timer)
  @timeit timer "DirectSolver - setup" begin
    asm = StaticAssembler(obj.domain)
    update_unknown_dofs!(obj.domain, asm)
    R = asm.residuals[obj.domain.dof.unknown_dofs]
    K = SparseArrays.sparse!(asm) |> Symmetric
    assumptions = OperatorAssumptions(false)
    prob = LinearProblem(K, R; assumptions=assumptions)
    linsolve = init(prob)
  end
  return DirectSolver(asm, linsolve, timer)
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
  @timeit timer(solver) "DirectSolver - gradient!" begin
    R = solver.assembler.residuals
    # annoying below TODO fix this
    update_field_unknowns!(p.U, obj.domain, Uu)
    R .= zero(eltype(R))
    domain_iterator!(R, obj.gradient, obj.domain, Uu, p)
    # set filled residual in linear solve init
    solver.linsolve.b = R[obj.domain.dof.unknown_dofs]
  end
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function hessian!(solver::DirectSolver, obj, Uu, p)
  @timeit timer(solver) "DirectSolver - hessian!" begin
    K = hessian!(solver.assembler, obj, Uu, p)
    solver.linsolve.A = K
  end
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function residual_norm(solver::DirectSolver, obj, Uu, p)
  # gradient!(solver, obj, Uu, p)
  return norm(solver.linsolve.b)
end

"""
$(TYPEDSIGNATURES)
"""
function solve!(ΔUu, solver::DirectSolver, obj, Uu, p)
  @timeit timer(solver) "DirectSolver - solve!" begin
    gradient!(solver, obj, Uu, p)
    hessian!(solver, obj, Uu, p)
    sol = LinearSolve.solve!(solver.linsolve)  
    ΔUu .= sol.u
  end
  return nothing
end

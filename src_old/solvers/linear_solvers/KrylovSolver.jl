mutable struct ObjectiveOperator{T}
  op::T
end

function ObjectiveOperator(obj::Objective, Uu, p)
  n = length(Uu)
  op = FunctionOperator(
    (du, u, p, t) -> hvp!(du, obj, Uu, p, u), zeros(n), zeros(n);
    p=p,
    t=0.0,
    isinplace=true
  )
  # return ObjectiveOperator(op)
  return op
end
# Base.size(A::ObjectiveOperator) = size(A.op)

function update(obj::Objective, Uu, p)
  n = length(Uu)
  # A.op = LinearOperator(
  #   Float64, n, n, true, false,
  #   (res, v, α, β) -> hvp!(res, obj, Uu, p, v)
  # )
  A = FunctionOperator(
    (du, u, p, t) -> hvp!(du, obj, Uu, p, u), zeros(n), zeros(n);
    p=p,
    t=p.t.current_time,
    isinplace=true
  )
  # return nothing
  return A
end

struct KrylovSolver{O, R, S, T} <: AbstractLinearSolver
  operator::O
  residual::R
  solver::S
  timer::T
end
timer(s::KrylovSolver) = s.timer

function KrylovSolver(obj::Objective, p, timer)
  update_unknown_dofs!(obj.domain)
  Uu = create_unknowns(obj.domain)
  R = create_fields(obj.domain)
  A = ObjectiveOperator(obj, Uu, p)
  # b = R[obj.domain.unknown_dofs]
  solver = CgSolver(length(Uu), length(Uu), typeof(Uu))
  # solver = GmresSolver(length(Uu), length(Uu), 20, typeof(Uu))
  return KrylovSolver(A, R, solver, timer)
end

function gradient!(solver::KrylovSolver, obj::Objective, Uu, p::ObjectiveParameters)
  @timeit timer(solver) "KrylovSolver - gradient!" begin
    update_field_unknowns!(p.U, obj.domain, Uu)
    solver.residual .= zero(eltype(solver.residual))
    domain_iterator!(solver.residual, obj.gradient, obj.domain, Uu, p)
  end
  return nothing
end

function solve!(ΔUu, solver::KrylovSolver, obj, Uu, p)
  @timeit timer(solver) "KrylovSolver - solve!" begin
    # update!(solver.operator, obj, Uu, p)
    A = update(obj, Uu, p)
    Pl = SciMLOperators.InvertedOperator(A)
    gradient!(solver.residual, obj, Uu, p)
    # A = solver.operator.op
    b = solver.residual[obj.domain.dof.unknown_dofs]
    # Pl = Preconditioners.AMGPreconditioner{RugeStuben}(A)
    # Pl = (A)
    # Pl = ldlt(A)
    # warm_start!(solver.solver, Uu)
    cg!(solver.solver, A, b, M=Pl, ldiv=true)
    # gmres!(solver.solver, A, b, M=Pl, ldiv=true)
  end
  ΔUu .= solver.solver.x
  return nothing
end

function residual_norm(solver::KrylovSolver, obj, Uu, p)
  return norm(solver.residual[obj.domain.dof.unknown_dofs])
end

abstract type AbstractLinearSolver end
abstract type AbstractNonlinearSolver end
function setup_solver_cache end

function Base.show(io::IO, solver::T) where T <: AbstractNonlinearSolver
  println(io, T.name.name)
  println(io, add_tabs(solver.linear_solver, 2))
end

# Direct linear solver implementation
struct DirectSolver{L} <: AbstractLinearSolver
  linsolve::L
end

function Base.similar(solver::DirectSolver)
  return DirectSolver(deepcopy(solver.linsolve))
end

function Base.show(io::IO, sol::DirectSolver)
  println(io, "DirectSolver:")
  println(io, "  Algorithm            = $(sol.linsolve.alg)")
  println(io, "  Assumptions          = $(sol.linsolve.assumptions)")
  println(io, "  Left preconditioner  = $(sol.linsolve.Pl)")
  println(io, "  Right preconditioner = $(sol.linsolve.Pr)")
end

function DirectSolver(domain::Domain)
  R = domain.solver_cache.assembler.residuals[domain.dof.unknown_dofs]
  K = SparseArrays.sparse!(domain.solver_cache) |> Symmetric
  assumptions = OperatorAssumptions(true)
  prob = LinearProblem(K, R; assumptions=assumptions)
  linsolve = init(prob)
  return DirectSolver{typeof(linsolve)}(linsolve)
end

function solve!(
  solver::DirectSolver, obj, domain, Uu
  # solver_alg
)
  gradient!(obj, domain, Uu)
  solver.linsolve.b = domain.solver_cache.assembler.residuals[domain.dof.unknown_dofs]
  solver.linsolve.A = hessian!(obj, domain, Uu) |> Symmetric
  # LinearSolve.solve!(solver.linsolve, solver_alg)
  # LinearSolve.solve!(solver.linsolve, CHOLMODFactorization())
  LinearSolve.solve!(solver.linsolve, CHOLMODFactorization())
  return nothing
end

function solve(solver::DirectSolver, obj, domain::Domain, Uu)
  solve!(solver, obj, domain, Uu)
  return solver.linsolve.u
end

# Matrix free linear solver implementation
struct MatrixFreeOperator{T, Field} <: LinearMaps.LinearMap{T}
  domain::Domain
  V::Field
  size::Dims{2}
end

function MatrixFreeOperator(domain::Domain)
  V = create_fields(domain)
  num_total_dofs = length(domain.dof.unknown_dofs)
  return MatrixFreeOperator{Float64, typeof(V)}(domain, V, (num_total_dofs, num_total_dofs))
end

Base.size(op::MatrixFreeOperator) = op.size

function LinearMaps._unsafe_mul!(Hv, H::MatrixFreeOperator, Vv::AbstractVector)
  # FiniteElementContainers.update_unknowns!(H.V, H.domain.dof, Vv)
  domain = H.domain
  stiffness_action!(domain, Vv)
  Hv .= domain.solver_cache.Hv[domain.dof.unknown_dofs]
  return nothing
end

struct MatrixFreeSolver{Op}
  op::Op
end

function MatrixFreeSolver(domain::Domain) 
  op = MatrixFreeOperator(domain)
  return MatrixFreeSolver(op)
end



# Newton Solver implementation
struct NewtonSolver{LinSolver} <: AbstractNonlinearSolver
  linear_solver::LinSolver
  max_iterations::Int
  absolute_tolerance::Float64
  relative_tolerance::Float64
end

function NewtonSolver(
  domain::Domain, linear_solver_type;
  max_iterations::Int = 1,
  absolute_tolerance::Float64 = 1e-8,
  relative_tolerance::Float64 = 1e-10
)
  linear_solver = linear_solver_type(domain)
  return NewtonSolver(linear_solver, max_iterations, absolute_tolerance, relative_tolerance)
end

function NewtonSolver(inputs::Dict{Symbol, Any})
  linear_solver_inputs = inputs[Symbol("linear solver")]
  linear_solver_types = keys(linear_solver_inputs) |> collect
  @assert length(linear_solver_types) == 1
  linear_solver_type = linear_solver_types[1]
  linear_solver = eval(linear_solver_type)(linear_solver_inputs[linear_solver_type])
  max_it = input_with_default(inputs, "maximum iterations", 10)
  abs_tol = input_with_default(inputs, "absolute tolerance", 1.0e-8)
  rel_tol = input_with_default(inputs, "relative tolerance", 1.0e-8)
  return NewtonSolver(linear_solver, max_it, abs_tol, rel_tol)
end

function Base.similar(solver::NewtonSolver)
  return NewtonSolver(
    similar(solver.linear_solver),
    solver.max_iterations,
    solver.absolute_tolerance,
    solver.relative_tolerance
  )
end

function logger(::NewtonSolver, n, norm_R, norm_R0, norm_U)
  @info @sprintf "  Iteration %5i: ||R|| = %1.6e    ||R/R0|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R (norm_R / norm_R0) norm_U
end

function solve(solver::NewtonSolver, obj, domain)
  Uu = create_unknowns(domain)
  solve!(Uu, solver, obj, domain)
  return Uu
end

function solve!(
  Uu, solver::NewtonSolver, obj, domain
  # linear_solver_alg
)
  @info "Newton Solver:"
  gradient!(obj, domain, Uu)
  R0 = norm(domain.solver_cache.assembler.residuals[domain.dof.unknown_dofs])
  for n in 1:10
    # solve!(solver.linear_solver, obj, domain, Uu, linear_solver_alg)
    solve!(solver.linear_solver, obj, domain, Uu)
    Uu .-= solver.linear_solver.linsolve.u
    norm_R = norm(solver.linear_solver.linsolve.b)
    norm_U = norm(solver.linear_solver.linsolve.u)
    logger(solver, n, norm_R, R0, norm_U)
    if norm_R < solver.absolute_tolerance || norm_U < solver.absolute_tolerance
      @info "Converged on absolute tolerance."
      break
    end
  end
  return nothing
end

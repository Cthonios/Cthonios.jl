struct LinearSolver{Solver, Preconditioner, Cache}
  solver::Solver
  precond::Preconditioner
  solver_cache::Cache
end

# TODO problably need some init methods

# TODO might need some constructors or something
function LinearSolver(
  input_settings::D, domain::QuasiStaticDomain
) where D <: Dict
  # make an initial dummy A and b
  Uu = create_unknowns(domain)
  # TODO need to set up custom solvers
  R = residual(domain, Uu)
  K = stiffness(domain, Uu)

  # setup solver
  if :type in keys(input_settings)
    if input_settings[:type] == "default"
      solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
    else
      solver = eval(Meta.parse(input_settings[:type]))()
    end
  else
    solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
  end

  # set up preconditioner
  if :preconditioner in keys(input_settings)
    if input_settings[:preconditioner] == "default"
      Pl = IdentityOperator(length(Uu))
    elseif input_settings[:preconditioner] == "cholesky"
      Pl = cholesky(K)
    else
      Pl = eval(Meta.parse(input_settings[:preconditioner]))(K)
    end
  else
    Pl = IdentityOperator(length(Uu))
  end

  prob   = LinearProblem(K, R)
  solver_cache = init(prob, solver; Pl=Pl) # TODO maybe we want a right preconditioner in some cases?

  return LinearSolver(solver, Pl, solver_cache)
end

function Base.show(io::IO, solver::LinearSolver)
  print(io, "LinearSolver\n", 
        "    Algorithm = $(typeof(solver.solver_cache.alg))")
end

solve(solver::LinearSolver) = LinearSolve.solve(solver.solver_cache)
solve!(solver::LinearSolver) = LinearSolve.solve!(solver.solver_cache)

# function update_residual!(
#   solver::LinearSolver, domain::QuasiStaticDomain, Uu::V
# ) where V <: AbstractVector
#   R = residual(domain, Uu)
#   solver.solver_cache.b = R
# end

# function update_stiffness!(
#   solver::LinearSolver, domain::QuasiStaticDomain, Uu::V
# ) where V <: AbstractVector
#   K = stiffness(domain, Uu)
#   solver.solver_cache.A = K
# end

# function update_residual_and_stiffness!(
#   solver::LinearSolver, domain::QuasiStaticDomain, Uu::V
# ) where V <: AbstractVector
#   update_residual!(solver, domain, Uu)
#   update_stiffness!(solver, domain, Uu)
# end

function update_residual!(
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}
  residual!(domain, Uu, domain.coords, U)
  solver.solver_cache.b .= @views domain.assembler.residuals[domain.dof.unknown_dofs]
end

function update_stiffness!(
  # K,
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}
  stiffness!(domain, Uu, domain.coords, U)
  solver.solver_cache.A = sparse(domain.assembler)
end



# function update!(solver::LinearSolver, domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
#   if solver.use_ad
#     @assert false "This should really happen in an extension"
#     # backend = AD.ForwardDiffBackend() # TODO customize AD backend
#     # # R = grad_energy_u(backend, domain, domain.coords, Uu)
#     # R = grad_energy_u(backend, domain, Uu, [], domain.time.current_time)
#     # hvp = hvp_energy_u(backend, domain, domain.coords, Uu)
#     # # f = (u, p, t) -> hvp(u)[1]
#     # f = FunctionWrapper{Vector{Float64}, Tuple{Vector{Float64}, Vector{Float64}, Float64}}((u, p, t) -> hvp(u)[1])
#     # K = FunctionOperator(f, zeros(length(Uu)), zeros(length(Uu)))
#     # solver.solver.A = K
#   else
#     R = residual(domain, Uu)
#     K = stiffness(domain, Uu)
#     solver.solver_cache.A = K
#   end
#   solver.solver_cache.b = R
# end


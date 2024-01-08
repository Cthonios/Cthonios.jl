struct LinearSolver{Solver, Preconditioner, Cache}
  solver::Solver
  precond::Preconditioner
  solver_cache::Cache
  use_ad::Bool
end

# TODO problably need some init methods

# TODO might need some constructors or something
function LinearSolver(input_settings::D, domain::QuasiStaticDomain) where D <: Dict
  # make an initial dummy A and b
  Uu = create_unknowns(domain)
  use_ad = input_settings["AD"]
  # TODO need to set up custom solvers
  if use_ad
    @info "Using AD to set up linear solver"
    @assert false "This should really happen in an extension"
    # backend = AD.ForwardDiffBackend() # TODO customize AD backend
    # R = grad_energy_u(backend, domain, Uu, [], domain.time.current_time)
    # # hvp = hvp_energy_u(backend, domain, domain.coords, Uu)
    # # hvp = hvp_energy_u(backend, domain, Uu)
    # hvp = hvp_energy_u(backend, domain, Uu, [], domain.time.current_time)
    # f = (u, p, t) -> hvp(u)[1]
    # # f = FunctionWrapper{Vector{Float64}, Tuple{Vector{Float64}, Vector{Float64}, Float64}}((u, p, t) -> hvp(u)[1])
    # K = FunctionOperator(f, zeros(length(Uu)), zeros(length(Uu)))
    # # Set up preconditioner here TODO
  else
    @info "Not using AD to set up linear solver, falling back to assembling K"
    R = residual(domain, Uu)
    K = stiffness(domain, Uu)
  end

  # setup solver
  if "type" in keys(input_settings)
    if input_settings["type"] == "default"
      solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
    else
      solver = eval(Meta.parse(input_settings["type"]))()
    end
  else
    solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
  end

  # set up preconditioner
  if "preconditioner" in keys(input_settings)
    if input_settings["preconditioner"] == "default"
      Pl = IdentityOperator(length(Uu))
    elseif input_settings["preconditioner"] == "cholesky"
      Pl = cholesky(K)
    else
      Pl = eval(Meta.parse(input_settings["preconditioner"]))(K)
    end
  else
    Pl = IdentityOperator(length(Uu))
  end

  prob   = LinearProblem(K, R)
  solver_cache = init(prob, solver; Pl=Pl) # TODO maybe we want a right preconditioner in some cases?

  return LinearSolver(solver, Pl, solver_cache, use_ad)
end

function Base.show(io::IO, solver::LinearSolver)
  print(io, "LinearSolver\n", 
        "    Algorithm = $(typeof(solver.solver_cache.alg))")
end

LinearSolve.solve(solver::LinearSolver) = solve(solver.solver_cache)

function update_residual!(solver::LinearSolver, domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
  R = residual(domain, Uu)
  solver.solver_cache.b = R
end

function update_stiffness!(solver::LinearSolver, domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
  K = stiffness(domain, Uu)
  solver.solver_cache.A = K
end

function update_residual_and_stiffness!(solver::LinearSolver, domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
  update_residual!(solver, domain, Uu)
  update_stiffness!(solver, domain, Uu)
end

function update_residual!(
  R::V1,
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}
  R .= zero(eltype(R))
  residual!(R, domain, Uu, U)
  solver.solver_cache.b = R[domain.dof.is_unknown]
end

function update_stiffness!(
  K,
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}
  K .= zero(eltype(K))
  stiffness!(K, domain, Uu, U)
  solver.solver_cache.A = K[domain.dof.is_unknown, domain.dof.is_unknown]
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

#########################################
# parsing
# TODO add a lot more settings, and make sure to make defaults
function read_linear_solvers(input_settings::D) where D <: Dict
  @assert "linear solvers" in keys(input_settings)
  solver_settings = input_settings["linear solvers"]
  solvers = Dict{String, Dict}()
  for solver_name in keys(solver_settings)
    @info "Reading linear solver $solver_name"
    settings = solver_settings[solver_name]

    @assert "type" in keys(settings)

    # will need to modify for things with preconditioners
    type = eval(Meta.parse(settings["type"]))

    if "preconditioner" in keys(settings)
      preconditioner = settings["preconditioner"]
    else
      preconditioner = Identity
    end

    if "AD" in keys(settings)
      @show "here"
      use_ad = settings["AD"]
    else
      use_ad = false
    end
    
    @info "  Type           = $type"
    @info "  Preconditioner = $preconditioner"
    @info "  AD             = $use_ad"

    solvers[solver_name] = Dict(
      "type"           => type,
      "preconditioner" => preconditioner,
      "AD"             => use_ad
    )
  end

  return solvers
end
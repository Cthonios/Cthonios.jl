struct LinearSolver{Solver}
  solver::Solver
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
    backend = AD.ForwardDiffBackend() # TODO customize AD backend
    R = grad_energy_u(backend, domain, domain.coords, Uu)
    hvp = hvp_energy_u(backend, domain, domain.coords, Uu)
    f = (u, p, t) -> hvp(u)[1]
    K = FunctionOperator(f, zeros(length(Uu)), zeros(length(Uu)))
    # Set up preconditioner here TODO
  else
    @info "Not using AD to set up linear solver, falling back to assembling K"
    R = residual(domain, Uu)
    K = stiffness(domain, Uu)
    # Set up preconditioner here TODO
  end
  prob   = LinearProblem(K, R) # TOOD set up preconditioner
  solver = init(prob)

  return LinearSolver(solver, use_ad)
end

function Base.show(io::IO, solver::LinearSolver)
  print(io, "LinearSolver\n", 
        "    Algorithm = $(typeof(solver.solver.alg))")
end

function update!(solver::LinearSolver, domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
  if solver.use_ad
    backend = AD.ForwardDiffBackend() # TODO customize AD backend
    R = grad_energy_u(backend, domain, domain.coords, Uu)
    hvp = hvp_energy_u(backend, domain, domain.coords, Uu)
    f = (u, p, t) -> hvp(u)[1]
    K = FunctionOperator(f, zeros(length(Uu)), zeros(length(Uu)))
  else
    R = residual(domain, Uu)
    K = stiffness(domain, Uu)
  end
  solver.solver.b = R
  solver.solver.A = K
end

function update_residual!(solver::LinearSolver, domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
  if solver.use_ad
    backend = AD.ForwardDiffBackend()
    R = grad_energy_u(backend, domain, domain.coords, Uu)
  else
    R = residual(domain, Uu)
  end
  solver.solver.b = R
end

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
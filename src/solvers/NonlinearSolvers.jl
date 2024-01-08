include("LinearSolvers.jl")

abstract type NonlinearSolverSettings end
abstract type NonlinearSolver{S, L, V} end

# for changing size after bc changes
function Base.resize!(solver::NonlinearSolver, domain::QuasiStaticDomain)
  resize!(solver.Uu, length(domain.dof.unknown_indices))
  resize!(solver.ΔUu, length(domain.dof.unknown_indices))
end
function update_increment!(solver::NonlinearSolver, ΔUu::V) where V <: AbstractVector
  solver.ΔUu .= ΔUu
end
function logger end # TO be defined for each solver
function solve! end # TO be defined for each solver

############
# Trust region solver
include("TrustRegionSolver.jl")

########################################################
# Newton solver
@kwdef struct NewtonSolverSettings <: NonlinearSolverSettings
  max_steps::Int = 10
end

Base.show(io::IO, settings::NewtonSolverSettings) = 
print(io, "NewtonSolverSettings\n", "    Max iterations = $(settings.max_steps)\n")

struct NewtonSolver{
  S <: NewtonSolverSettings, 
  L <: LinearSolver, 
  V <: AbstractVector
} <: NonlinearSolver{S, L, V}
  settings::S
  linear_solver::L
  Uu::V
  ΔUu::V
end

# TODO maybe add preconditioner below?
function NewtonSolver(input_settings::D, domain::QuasiStaticDomain) where D <: Dict
  settings      = NewtonSolverSettings() # TODO add non-defaults
  linear_solver = LinearSolver(input_settings["linear solver"], domain)
  Uu            = create_unknowns(domain)
  ΔUu           = create_unknowns(domain)
  return NewtonSolver(settings, linear_solver, Uu, ΔUu)
end

# function NewtonSolver()

function Base.show(io::IO, solver::NewtonSolver)
  print(io, "NewtonSolver\n", 
        "  Settings      = $(solver.settings)\n",
        "  Linear solver = $(solver.linear_solver)\n")
end

function logger(::NewtonSolver, n, norm_R, norm_U)
  @info @sprintf "Iteration %5i: ||R|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R norm_U
end

function solve!(solver::NewtonSolver, domain::QuasiStaticDomain)
  # unpack cached arrays from solver
  Uu, ΔUu = solver.Uu, solver.ΔUu
  U = create_fields(domain)
  R = create_fields(domain)
  # R = domain.assembler.R
  K = domain.assembler.K
  # update linear solver stiffness
  update_bcs!(U, domain.coords, domain.time.current_time, domain.bcs)
  # update_stiffness!(solver.linear_solver, domain, Uu)
  update_stiffness!(K, solver.linear_solver, domain, Uu, U)

  for n in 1:solver.settings.max_steps
    # update_residual!(solver.linear_solver, domain, Uu)
    update_residual!(R, solver.linear_solver, domain, Uu, U)
    sol = solve(solver.linear_solver)
    update_increment!(solver, sol.u)

    # # # TODO above should use linear solver in solver
    @. Uu    = Uu - ΔUu
    # norm_R   = norm(R)
    norm_R   = norm(solver.linear_solver.solver_cache.b)
    norm_ΔUu = norm(ΔUu)
    
    logger(solver, n, norm_R, norm_ΔUu)

    if norm_R <= 1e-12 || norm_ΔUu <= 1e-12
      break
    end
  end
end

################################
# parsing
function read_nonlinear_solvers(input_settings::D) where D <: Dict
  @assert "nonlinear solvers" in keys(input_settings)
  solver_settings = input_settings["nonlinear solvers"]
  solvers = Dict{String, Dict}()
  for solver_name in keys(solver_settings)
    settings = solver_settings[solver_name]
    @assert "type" in keys(settings)
    @assert "linear solver" in keys(settings)
    type = settings["type"]
    linear_solver = settings["linear solver"]
    @assert linear_solver in keys(input_settings["linear solvers"])
    solvers[solver_name] = Dict(
      "type" => type,
      "linear solver" => input_settings["linear solvers"][linear_solver]
    )
  end 
  return solvers
end

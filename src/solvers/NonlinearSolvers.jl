include("LinearSolvers.jl")

abstract type NonlinearSolverSettings end
abstract type NonlinearSolver{S, L, V} end

# for changing size after bc changes
function Base.resize!(solver::NonlinearSolver, domain::QuasiStaticDomain)
  resize!(solver.Uu, length(domain.dof.unknown_indices))
  resize!(solver.ΔUu, length(domain.dof.unknown_indices))
end
function logger end # TO be defined for each solver
function solve! end # TO be defined for each solver

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
  
  # update linear solver
  update!(solver.linear_solver, domain, Uu)

  for n in 1:solver.settings.max_steps
    # R = residual(domain, Uu)
    # K = stiffness(domain, Uu)
    # ΔUu .= K \ R
    # @show size(Uu)
    # @show size(solver.linear_solver.solver.b)
    # @show size(solver.linear_solver.solver.A)
    update_residual!(solver.linear_solver, domain, Uu)
    sol = solve(solver.linear_solver.solver)
    # LinearSolve.solve!(solver.linear_solver.solver)
    # return solver.linear_solver.solver
    ΔUu .= sol.u
    # ΔUu .= solver.linear_solver.solver.u
    # # @show ΔUu
    # # @show sol
    # # return sol

    # # # TODO above should use linear solver in solver
    @. Uu    = Uu - ΔUu
    # norm_R   = norm(R)
    norm_R   = norm(solver.linear_solver.solver.b)
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



# function solve!(domain::Domain, solver::NewtonSolver)
#   # op = QuasiStaticDomainOperator{Float64}(domain)
#   dof       = domain.dof
#   assembler = domain.assembler

#   # extra assembly for now since it check pos def at beginning
#   # update_fields!(domain)
#   # assemble!(domain)

#   for n in 1:solver.settings.max_steps
#     update_fields!(domain)
#     assemble!(domain)

#     K = assembler.K[dof.is_unknown, dof.is_unknown]
#     R = assembler.R[dof.is_unknown]
#     # UpdatePreconditioner!(P, K, 2)

#     # P = CholeskyPreconditioner(K, 2)
#     # P = AMGPreconditioner{RugeStuben}(K)
#     # P = AMGPreconditioner{SmoothedAggregation}(K)
#     # TODO genrealize solver
#     # cg!(solver.ΔUu, K, R; Pl=P)
#     # gmres!(solver.ΔUu, K, R; Pl=P)
#     # update unknowns

#     solver.ΔUu .= K \ R

#     domain.Uu .= domain.Uu .- solver.ΔUu

#     norm_R = norm(domain.assembler.R[domain.dof.is_unknown])
#     norm_U = norm(solver.ΔUu)

#     # TODO redo
#     @info @sprintf "Iteration %5i: ||R|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R norm_U

#     # TODO make solver tolerance a setting
#     if norm_R < 1e-12 || norm_U < 1e-12
#       break
#     end
#   end

# end

# include("TrustRegionSolver.jl")
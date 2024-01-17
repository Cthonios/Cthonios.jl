@kwdef struct NewtonSolverSettings <: NonlinearSolverSettings
  max_steps::Int = 20
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
  linear_solver = LinearSolver(input_settings[Symbol("linear solver")], domain)
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

function solve!(
  solver::NewtonSolver, domain::QuasiStaticDomain,
  common::CthoniosCommon
)
  # unpack cached arrays from solver
  Uu, ΔUu = solver.Uu, solver.ΔUu
  U = create_fields(domain)

  @timeit timer(common) "Update BCs" update_bcs!(U, domain, domain.coords)
  # @timeit timer(common) "Stiffness" update_stiffness!(solver.linear_solver, domain, Uu, U)

  for n in 1:solver.settings.max_steps
    # update_residual!(solver.linear_solver, domain, Uu)
    @timeit timer(common) "Residual" update_residual!(solver.linear_solver, domain, Uu, U)
    @timeit timer(common) "Stiffness" update_stiffness!(solver.linear_solver, domain, Uu, U)

    # @timeit timer(common) "Linear solve" sol = solve(solver.linear_solver)
    # update_increment!(solver, sol.u)


    # @timeit timer(common) "Linear solve" solve!(solver.linear_solver)
    # update_increment!(solver, solver.linear_solver.solver_cache.u)
    @timeit timer(common) "Linear solve" solve!(solver.linear_solver)
    update_increment!(solver, solver.linear_solver.unknowns)


    # # # TODO above should use linear solver in solver
    @. Uu    = Uu - ΔUu
    # norm_R   = norm(R)
    # norm_R   = norm(solver.linear_solver.solver_cache.b)
    norm_R   = norm(solver.linear_solver.residual)
    norm_ΔUu = norm(ΔUu)
    
    logger(solver, n, norm_R, norm_ΔUu)

    if norm_R <= 1e-12 || norm_ΔUu <= 1e-12
      break
    end
  end
end

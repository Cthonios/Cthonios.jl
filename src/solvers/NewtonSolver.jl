@kwdef struct NewtonSolverSettings <: NonlinearSolverSettings
  max_steps::Int = 20
  relative_tolerance::Float64 = 1.0e-8
  absolute_tolerance::Float64 = 1.0e-8
  update_stiffness_each_iteration::Bool = false
end

function NewtonSolverSettings(input_settings::D) where D <: Dict
  if Symbol("absolute tolerance") in keys(input_settings)
    absolute_tolerance = input_settings[Symbol("absolute tolerance")]
  else
    absolute_tolerance = 1.0e-8
  end

  if Symbol("max steps") in keys(input_settings)
    max_steps = input_settings[Symbol("max steps")]
  else
    max_steps = 20
  end

  if Symbol("relative tolerance") in keys(input_settings)
    relative_tolerance = input_settings[Symbol("relative tolerance")]
  else
    relative_tolerance = 1.0e-8
  end

  if Symbol("update stiffness each iteration") in keys(input_settings)
    update_stiffness_each_iteration = input_settings[Symbol("update stiffness each iteration")]
  else
    update_stiffness_each_iteration = false
  end

  return NewtonSolverSettings(
    max_steps, 
    relative_tolerance, absolute_tolerance, 
    update_stiffness_each_iteration
  )
end

function Base.show(io::IO, settings::NewtonSolverSettings)
  print(io, "NewtonSolverSettings\n", 
            "  Max iterations                  = $(settings.max_steps)\n",
            "  Relative tolerance              = $(settings.relative_tolerance)\n",
            "  Absolute tolerance              = $(settings.absolute_tolerance)\n",
            "  Update stiffness each iteration = $(settings.update_stiffness_each_iteration)\n"
  )
end

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
  settings      = NewtonSolverSettings(input_settings) # TODO add non-defaults
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

function logger(::NewtonSolver, n, norm_R, norm_R0, norm_U)
  @info @sprintf "Iteration %5i: ||R|| = %1.6e    ||R/R0|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R (norm_R / norm_R0) norm_U
end

function solve!(
  solver::NewtonSolver, domain::QuasiStaticDomain,
  common::CthoniosCommon
)
  # unpack cached arrays from solver
  Uu, ΔUu = solver.Uu, solver.ΔUu
  U = create_fields(domain)

  @timeit timer(common) "Update BCs" update_bcs!(U, domain, domain.coords)
  @timeit timer(common) "Stiffness" update_stiffness!(solver.linear_solver, domain, Uu, U, common)

  norm_R0 = 0.0
  for n in 1:solver.settings.max_steps
    # update_residual!(solver.linear_solver, domain, Uu)
    @timeit timer(common) "Residual" update_residual!(solver.linear_solver, domain, Uu, U)

    if solver.settings.update_stiffness_each_iteration && n > 2
      @timeit timer(common) "Stiffness" update_stiffness!(solver.linear_solver, domain, Uu, U, common)
    end

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
    
    if n == 1
      norm_R0 = norm_R
    end

    logger(solver, n, norm_R, norm_R0, norm_ΔUu)

    # if norm_R <= 1e-12 || norm_ΔUu <= 1e-12
    if norm_R   <= solver.settings.absolute_tolerance ||
       norm_ΔUu <= solver.settings.absolute_tolerance
      break
    end

    if (norm_R / norm_R0) <= solver.settings.relative_tolerance
      break
    end
  end
end

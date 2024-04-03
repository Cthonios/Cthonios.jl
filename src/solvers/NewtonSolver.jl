@kwdef struct NewtonSolverSettings <: NonlinearSolverSettings
  max_steps::Int = 20
  relative_tolerance::Float64 = 1.0e-8
  absolute_tolerance::Float64 = 1.0e-8
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

  return NewtonSolverSettings(
    max_steps, 
    relative_tolerance, absolute_tolerance
  )
end

function Base.show(io::IO, settings::NewtonSolverSettings)
  print(io, "NewtonSolverSettings\n", 
            "        Max iterations                  = $(settings.max_steps)\n",
            "        Relative tolerance              = $(settings.relative_tolerance)\n",
            "        Absolute tolerance              = $(settings.absolute_tolerance)\n"
  )
end

struct NewtonSolver{
  S <: NewtonSolverSettings, 
  L <: AbstractLinearSolver
} <: NonlinearSolver{S, L}
  settings::S
  linear_solver::L
  use_warm_start::Bool
end

# TODO maybe add preconditioner below?
function NewtonSolver(input_settings::D, domain::QuasiStaticDomain, backend) where D <: Dict
  settings      = NewtonSolverSettings(input_settings) # TODO add non-defaults
  linear_solver = setup_linear_solver(input_settings[Symbol("linear solver")], domain, backend)

  if Symbol("warm start") in keys(input_settings)
    parsed_string = input_settings[Symbol("warm start")]
    if parsed_string == "on"
      use_warm_start = true
    elseif parsed_string == "off"
      use_warm_start = false
    else
      @assert false "Bad value for warm start"
    end
  else
    use_warm_start = false
  end
  return NewtonSolver(settings, linear_solver, use_warm_start)
end

function Base.show(io::IO, solver::NewtonSolver)
  print(io, "\n    NewtonSolver\n", 
        "      $(solver.settings)\n",
        "    $(solver.linear_solver)\n")
end

function logger(::NewtonSolver, n, norm_R, norm_R0, norm_U)
  @info @sprintf "Iteration %5i: ||R|| = %1.6e    ||R/R0|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R (norm_R / norm_R0) norm_U
end

function solve!(
  solver::NewtonSolver, domain::QuasiStaticDomain,
  common::CthoniosCommon
)
  # unpack cached arrays from solver and domain
  @unpack X, Uu, ΔUu, U = domain.domain_cache

  @timeit timer(common) "Update BCs" begin 
    update_bcs!(U, domain, X)
  end

  norm_R0 = 0.0
  for n in 1:solver.settings.max_steps
    @timeit timer(common) "Residual and stiffness" begin 
      internal_force_and_stiffness!(solver, domain, domain.domain_cache, Uu, backend(common))
    end

    @timeit timer(common) "Linear solve" begin 
      solve!(ΔUu, solver.linear_solver, domain, common)
    end

    # TODO above should use linear solver in solver
    @. Uu    = Uu - ΔUu
    # norm_R   = @views norm(domain.domain_cache.f[domain.dof.unknown_dofs])
    norm_R   = @views norm(solver.linear_solver.assembler.residuals[domain.dof.unknown_dofs])
    norm_ΔUu = norm(ΔUu)
    
    # save first residual norm for calculating relative residual
    if n == 1
      norm_R0 = norm_R
    end

    logger(solver, n, norm_R, norm_R0, norm_ΔUu)

    # check convergence on absolute tolerance
    if norm_R   <= solver.settings.absolute_tolerance ||
       norm_ΔUu <= solver.settings.absolute_tolerance
      return nothing
    end

    # check convergence on relative tolerance
    if (norm_R / norm_R0) <= solver.settings.relative_tolerance
      return nothing
    end
  end

  max_newton_steps_exception()
end

struct MaxNewtonStepsException <: Exception
end

function Base.show(io::IO, ::MaxNewtonStepsException)
  println(io, "Maximum number of newton steps reached!")
end

function max_newton_steps_exception()
  throw(MaxNewtonStepsException())
end
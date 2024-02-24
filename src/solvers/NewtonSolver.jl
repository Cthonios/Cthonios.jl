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
            "        Max iterations                  = $(settings.max_steps)\n",
            "        Relative tolerance              = $(settings.relative_tolerance)\n",
            "        Absolute tolerance              = $(settings.absolute_tolerance)\n",
            "        Update stiffness each iteration = $(settings.update_stiffness_each_iteration)\n"
  )
end

struct NewtonSolver{
  S <: NewtonSolverSettings, 
  L <: AbstractLinearSolver, 
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
  linear_solver = setup_linear_solver(input_settings[Symbol("linear solver")], domain)
  Uu            = create_unknowns(domain)
  ΔUu           = create_unknowns(domain)
  return NewtonSolver(settings, linear_solver, Uu, ΔUu)
end

function Base.show(io::IO, solver::NewtonSolver)
  print(io, "\n    NewtonSolver\n", 
        "      $(solver.settings)\n",
        "    $(solver.linear_solver)\n")
end

function logger(::NewtonSolver, n, norm_R, norm_R0, norm_U)
  @info @sprintf "Iteration %5i: ||R|| = %1.6e    ||R/R0|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R (norm_R / norm_R0) norm_U
end

function update_unknown_dofs!(solver::NewtonSolver, d::QuasiStaticDomain)
  # update the dofs
  FiniteElementContainers.update_unknown_dofs!(d.dof, d.bc_dofs)
  FiniteElementContainers.update_unknown_dofs!(solver.linear_solver.assembler, d.dof, map(x -> x.fspace, d.sections), d.bc_dofs)

  # now resize the caches
  resize!(solver.Uu, length(d.dof.unknown_dofs))
  resize!(solver.ΔUu, length(d.dof.unknown_dofs))
end

function solve!(
  solver::NewtonSolver, domain::QuasiStaticDomain,
  common::CthoniosCommon
)
  # unpack cached arrays from solver and domain
  Uu, ΔUu = solver.Uu, solver.ΔUu
  U = domain.domain_cache.U

  @timeit timer(common) "Update BCs" update_bcs!(U, domain, domain.coords)

  norm_R0 = 0.0
  for n in 1:solver.settings.max_steps
    @timeit timer(common) "Residual and stiffness" begin 
      update_residual_and_stiffness!(solver.linear_solver, domain, Uu, U)
    end

    @timeit timer(common) "Linear solve" solve!(ΔUu, solver.linear_solver, domain, common)

    # TODO above should use linear solver in solver
    @. Uu    = Uu - ΔUu
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
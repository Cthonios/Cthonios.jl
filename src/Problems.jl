abstract type AbstractProblem <: AbstractCthoniosType end


struct ForwardProblem{Domain, Solver}
  domain::Domain
  solver::Solver
end

function ForwardProblem(input_settings::D) where D <: Dict
  domain            = QuasiStaticDomain(input_settings)
  # update bcs once here, TODO this will fail when we start disabling bcs
  # we need to do this for the solver, and also need to probably reinit
  # thie linear solver when we change bcs
  update_unknown_ids!(domain)

  nonlinear_solvers = read_nonlinear_solvers(input_settings)
  @assert "problem" in keys(input_settings)
  @assert "solver" in keys(input_settings["problem"])
  solver_settings   = nonlinear_solvers[input_settings["problem"]["solver"]]
  solver_type       = Meta.parse(solver_settings["type"])
  solver            = eval(solver_type)(solver_settings, domain)
  return ForwardProblem(domain, solver)
end

ForwardProblem(f::String) = ForwardProblem(YAML.load_file(f))

function Base.show(io::IO, problem::ForwardProblem)
  println(io, "ForwardProblem")
  println(io, "  Domain", problem.domain)
  println(io, "  Solver", problem.solver)
end

function solve!(problem::ForwardProblem)
  domain, solver = problem.domain, problem.solver
  reset!(domain.time)

  # TODO need a way to update bcs properly, this will change a lot
  # update_unknown_ids!(domain)
  # resize!(solver, domain)

  post_process_load_step!(domain, solver)

  while domain.time.current_time <= domain.time.end_time
    step!(domain.time)

    # TODO need a way to update bcs properly, this will change a lot
    # update_unknown_ids!(domain)
    # resize!(solver, domain)

    @info "$(repeat('=', 64))"
    @info "= Load step $(domain.time.current_time_step - 1)"
    @info "= Time      $(domain.time.current_time)"
    @info "$(repeat('=', 64))"

    solve!(solver, domain)

    # post-processing
    post_process_load_step!(domain, solver) 
  end

  close(domain.post_processor)
end

function post_process_load_step!(domain::QuasiStaticDomain, solver::NonlinearSolver)
  # update_bcs!(domain.post_processor.scratch_U, domain)
  update_bcs!(domain.post_processor.scratch_U, domain.coords, domain.time.current_time, domain.bcs)
  update_fields!(domain.post_processor.scratch_U, domain, solver.Uu)

  # write displacements
  write_time(domain.post_processor, domain.time.current_time_step, domain.time.current_time)
  write_values(domain.post_processor, NodalVariable, domain.time.current_time_step, "displ_x", 
               domain.post_processor.scratch_U[1, :])
  write_values(domain.post_processor, NodalVariable, domain.time.current_time_step, "displ_y", 
               domain.post_processor.scratch_U[2, :])

  # TODO fill out more output
end
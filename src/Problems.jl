abstract type AbstractProblem <: AbstractCthoniosType end


struct ForwardProblem{Domain, Solver, Post}
  domain::Domain
  solver::Solver
  post_processor::Post
end

function get_domain_input_settings(inputs::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any} 
  return inputs[:domain]
end

function get_mesh_file_name(inputs::D)::String where D <: Dict{Symbol, Any}
  return inputs[:domain][:mesh][Symbol("file name")]
end

function get_output_file_name(inputs::D)::String where D <: Dict{Symbol, Any}
  return inputs[:results][Symbol("output file name")]
end

function get_solver_input_settings(inputs::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return inputs[:solver]
end

function ForwardProblem(input_settings::D, common::CthoniosCommon) where D <: Dict

  @timeit timer(common) "Setup" begin
    @timeit timer(common) "Domain" domain = QuasiStaticDomain(get_domain_input_settings(input_settings))
    # update bcs once here, TODO this will fail when we start disabling bcs
    # we need to do this for the solver, and also need to probably reinit
    # thie linear solver when we change bcs
    # @timeit timer(common) "Update unknown dofs" update_unknown_dofs!(domain)

    # set up solver
    @timeit timer(common) "Solver" solver = setup_nonlinear_solver(get_solver_input_settings(input_settings), domain)
    @timeit timer(common) "Update unknown dofs" update_unknown_dofs!(solver, domain)
    @timeit timer(common) "Postprocessor" post_processor = PostProcessor(
      get_mesh_file_name(input_settings), get_output_file_name(input_settings),
      domain.dof, size(domain.domain_cache.X, 1)
    )
  end
  return ForwardProblem(domain, solver, post_processor)
end

function Base.show(io::IO, problem::ForwardProblem)
  println(io, "ForwardProblem")
  println(io, "  Domain", problem.domain)
  println(io, "\n  Solver\n", problem.solver)
end

function solve!(problem::ForwardProblem, common::CthoniosCommon)
  domain, solver = problem.domain, problem.solver
  reset!(domain.time)

  # TODO need a way to update bcs properly, this will change a lot
  # update_unknown_ids!(domain)
  # resize!(solver, domain)

  @timeit timer(common) "Results output" post_process_load_step!(problem)

  while domain.time.current_time <= domain.time.end_time
    step!(domain.time)

    # TODO need a way to update bcs properly, this will change a lot
    # update_unknown_ids!(domain)
    # resize!(solver, domain)

    # @info "$(repeat('=', 64))"
    @info "$(repeat('=', 96))"
    @info "= Load step $(domain.time.current_time_step - 1)"
    @info "= Time      $(domain.time.current_time)"
    # @info "$(repeat('=', 64))"
    @info "$(repeat('=', 96))"

    @timeit timer(common) "Nonlinear solver" solve!(solver, domain, common)

    # post-processing
    @timeit timer(common) "Results output" post_process_load_step!(problem) 

  end

  close(problem.post_processor)
end

function post_process_load_step!(problem::ForwardProblem)
  # update_bcs!(domain.post_processor.scratch_U, domain)
  domain = problem.domain
  solver = problem.solver
  pp = problem.post_processor
  update_fields!(pp.scratch_U, domain, domain.domain_cache.X, solver.Uu)

  # write displacements
  write_time(pp, domain.time.current_time_step, domain.time.current_time)
  write_values(pp, NodalVariable, domain.time.current_time_step, "displ_x", pp.scratch_U[1, :])
  write_values(pp, NodalVariable, domain.time.current_time_step, "displ_y", pp.scratch_U[2, :])

  if size(problem.domain.domain_cache.X, 1) == 3
    write_values(pp, NodalVariable, domain.time.current_time_step, "displ_z", pp.scratch_U[3, :])
  end

  # TODO fill out more output
end
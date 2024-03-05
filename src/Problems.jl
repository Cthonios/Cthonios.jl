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

function get_output_nodal_fields(inputs::D)::Vector{String} where D <: Dict{Symbol, Any}
  return inputs[:results][Symbol("nodal fields")]
end

function get_output_element_fields(inputs::D)::Vector{String} where D <: Dict{Symbol, Any}
  if Symbol("element fields") in keys(inputs[:results])
    return inputs[:results][Symbol("element fields")]
  else
    return Vector{String}(undef, 0)
  end
end

function get_output_quadrature_fields(inputs::D)::Vector{String} where D <: Dict{Symbol, Any}
  if Symbol("quadrature fields") in keys(inputs[:results])
    return inputs[:results][Symbol("quadrature fields")]
  else
    return Vector{String}(undef, 0)
  end
end

function get_solver_input_settings(inputs::D)::Dict{Symbol, Any} where D <: Dict{Symbol, Any}
  return inputs[:solver]
end

function get_max_properties(domain)
  return maximum(map(x -> ConstitutiveModels.num_properties(x.model), domain.sections))
end

function get_max_state_variables(domain)
  return maximum(map(x -> ConstitutiveModels.num_state_vars(x.model), domain.sections))
end

function get_max_q_points(domain)
  return maximum(map(x -> FiniteElementContainers.num_q_points(x), domain.sections))
end

function ForwardProblem(input_settings::D, common::CthoniosCommon) where D <: Dict

  @timeit timer(common) "Setup" begin
    @timeit timer(common) "Domain" domain = QuasiStaticDomain(get_domain_input_settings(input_settings))
    # update bcs once here, TODO this will fail when we start disabling bcs
    # we need to do this for the solver, and also need to probably reinit
    # thie linear solver when we change bcs
    # @timeit timer(common) "Update unknown dofs" update_unknown_dofs!(domain)

    # set up solver
    @timeit timer(common) "Solver" solver = setup_nonlinear_solver(get_solver_input_settings(input_settings), domain, backend(common))
    @timeit timer(common) "Update unknown dofs" update_unknown_dofs!(solver, domain)
    @timeit timer(common) "Postprocessor" post_processor = PostProcessor(
      get_mesh_file_name(input_settings), 
      get_output_file_name(input_settings),
      get_output_nodal_fields(input_settings),
      get_output_element_fields(input_settings),
      get_output_quadrature_fields(input_settings),
      size(domain.domain_cache.X, 1),
      get_max_properties(domain), 
      get_max_state_variables(domain),
      get_max_q_points(domain)
    )
  end
  return ForwardProblem(domain, solver, post_processor)
end

function Base.show(io::IO, problem::ForwardProblem)
  println(io, "ForwardProblem")
  println(io, "  Domain", problem.domain)
  println(io, "\n  Solver\n", problem.solver)
  println(io, "\n  Results\n", problem.post_processor)
end

function solve!(problem::ForwardProblem, common::CthoniosCommon)
  domain, solver = problem.domain, problem.solver
  reset!(domain.domain_cache.time)

  @timeit timer(common) "Results output" post_process_load_step!(problem, common)

  while domain.domain_cache.time.current_time <= domain.domain_cache.time.end_time
    step!(domain.domain_cache.time)

    @info "$(repeat('=', 96))"
    @info "= Load step $(domain.domain_cache.time.current_time_step - 1)"
    @info "= Time      $(domain.domain_cache.time.current_time)"
    @info "$(repeat('=', 96))"

    @timeit timer(common) "Nonlinear solver" solve!(solver, domain, common)

    # copy new state variables to old array
    @timeit timer(common) "State variable update" begin
      domain.domain_cache.state_old .= domain.domain_cache.state_new
    end
    # post-processing
    @timeit timer(common) "Results output" post_process_load_step!(problem, common) 

  end

  close(problem.post_processor)
end

function post_process_load_step!(problem::ForwardProblem, common::CthoniosCommon)
  # update_bcs!(domain.post_processor.scratch_U, domain)
  domain = problem.domain
  cache = problem.domain.domain_cache
  solver = problem.solver
  pp = problem.post_processor
  time_step = domain.domain_cache.time.current_time_step
  # update_fields!(pp.scratch_U, domain, domain.domain_cache.X, solver.Uu)

  # write displacements
  write_time(pp, domain.domain_cache.time.current_time_step, domain.domain_cache.time.current_time)

  # Nodal variables first
  if "displ_x" in keys(pp.out_file.nodal_var_name_dict)
    write_values(pp, NodalVariable, time_step, "displ_x", cache.U[1, :])
    write_values(pp, NodalVariable, time_step, "displ_y", cache.U[2, :])

    if size(problem.domain.domain_cache.X, 1) == 3
      write_values(pp, NodalVariable, time_step, "displ_z", cache.U[3, :])
    end
  end

  if "internal_force_x" in keys(pp.out_file.nodal_var_name_dict)
    # need to run the routine here
    internal_force!(solver, domain, solver.Uu, backend(common))

    write_values(pp, NodalVariable, time_step, "internal_force_x", cache.f[1, :])
    write_values(pp, NodalVariable, time_step, "internal_force_y", cache.f[2, :])

    if size(problem.domain.domain_cache.X, 1) == 3
      write_values(pp, NodalVariable, time_step, "internal_force_z", cache.f[3, :])
    end
  end

  # Now element variables
  for (name, section) in pairs(domain.sections)
    if "properties_001" in keys(pp.out_file.element_var_name_dict)
      for n in axes(domain.domain_cache.props[name], 1)
        temp = @views cache.props[name]
        write_values(pp, ElementVariable, time_step, section.block_id, "properties_$(lpad(n, 3, '0'))", temp[n, :])
      end
    end
  end

  # TODO fill out more output
end

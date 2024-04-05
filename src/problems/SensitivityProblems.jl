abstract type AbstractSensitityProblem end

struct EnergySensitivityProblem{
  Domain,
  Solver,
  Post
} <: AbstractSensitityProblem
  domain::Domain
  solver::Solver
  post_processor::Post
end

function EnergySensitivityProblem(prob::ForwardProblem, input_settings, common)
  domain = prob.domain
  solver = prob.solver
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
  
  return EnergySensitivityProblem(domain, solver, post_processor)
end

function solve!(prob::EnergySensitivityProblem, common::CthoniosCommon)
  solver = prob.domain
  domain = prob.domain
  cache  = prob.domain.domain_cache
  tangent_cache = similar(cache)
  seed!(Reverse, tangent_cache)
  tangent_cache.Π .= 1.0
  Uu = cache.Uu
  dUu = similar(Uu)

  autodiff(
    Reverse, strain_energy!,
    Const(solver),
    Const(domain),
    Duplicated(cache, tangent_cache),
    Duplicated(Uu, dUu),
    Const(backend(common))
  )

  post_process!(prob, tangent_cache)
  close(prob.post_processor)
end

function post_process!(prob::EnergySensitivityProblem, tangent_cache)
  time_step = 1
  domain = prob.domain
  pp = prob.post_processor
  U  = prob.domain.domain_cache.U
  write_time(pp, 1, 0.0)

  if "dcoordinates_x" in keys(pp.out_file.nodal_var_name_dict)
    write_values(pp, NodalVariable, time_step, "dcoordinates_x", tangent_cache.X[1, :])
    write_values(pp, NodalVariable, time_step, "dcoordinates_y", tangent_cache.X[2, :])

    if size(prob.domain.domain_cache.X, 1) == 3
      write_values(pp, NodalVariable, time_step, "dcoordinates_z", tangent_cache.X[3, :])
    end
  end

  if "displ_x" in keys(pp.out_file.nodal_var_name_dict)
    write_values(pp, NodalVariable, time_step, "displ_x", U[1, :])
    write_values(pp, NodalVariable, time_step, "displ_y", U[2, :])

    if size(prob.domain.domain_cache.X, 1) == 3
      write_values(pp, NodalVariable, time_step, "displ_z", U[3, :])
    end
  end

  for (name, section) in pairs(domain.sections)
    if "dproperties_001" in keys(pp.out_file.element_var_name_dict)
      for n in axes(domain.domain_cache.props[name], 1)
        temp = @views tangent_cache.props[name]
        write_values(pp, ElementVariable, time_step, section.block_id, "dproperties_$(lpad(n, 3, '0'))", temp[n, :])
      end
    end
  end
end
abstract type AbstractProblem <: AbstractCthoniosType end


struct ForwardProblem{Domain, Solver}
  domain::Domain
  solver::Solver
end

function ForwardProblem(input_settings::D) where D <: Dict
  @assert "solver" in keys(input_settings)
  @assert "type"   in keys(input_settings["solver"])
  type = input_settings["solver"]["type"]
  # TODO generalize
  domain = QuasiStaticDomain(input_settings)
  # solver = NewtonSolver(input_settings, domain) 
  # solver = TrustRegionSolver(input_sett)
  solver = eval(Meta.parse(type))(input_settings, domain)
  return ForwardProblem(domain, solver)
end

Base.show(io::IO, ::ForwardProblem) = println(io, "Forward Problem - Fill this out")

function solve!(problem::ForwardProblem)
  domain, solver = problem.domain, problem.solver
  reset!(domain.time)
  domain.U  .= 0.0
  domain.Uu .= 0.0

  write_time(domain.post_processor, domain.time.current_time_step, 0.0)
  write_values(domain.post_processor, NodalVariable, domain.time.current_time_step, "displ_x", domain.U[1, :])
  write_values(domain.post_processor, NodalVariable, domain.time.current_time_step, "displ_y", domain.U[2, :])

  while domain.time.current_time <= domain.time.end_time
    # time update
    step!(domain.time)

    @info "$(repeat('=', 64))"
    @info "= Load step $(domain.time.current_time_step - 1)"
    @info "= Time      $(domain.time.current_time)"
    @info "$(repeat('=', 64))"
    # update bcs TODO wrap below in method
    # so we don't accidently not resize in future
    update_unknown_ids!(domain)
    resize!(domain.Uu, sum(domain.dof.is_unknown))
    resize!(solver.ΔUu, sum(domain.dof.is_unknown))
    solver.ΔUu .= 0.0 # Unitful issue
    update_bcs!(domain)
    update_fields!(domain)

    # solve
    solve!(domain, solver)
    update_fields!(domain)
    write_time(domain.post_processor, domain.time.current_time_step, domain.time.current_time)
    write_values(domain.post_processor, NodalVariable, domain.time.current_time_step, "displ_x", domain.U[1, :])
    write_values(domain.post_processor, NodalVariable, domain.time.current_time_step, "displ_y", domain.U[2, :])

    # step!(domain.time)
  end

  close(domain.post_processor)
end

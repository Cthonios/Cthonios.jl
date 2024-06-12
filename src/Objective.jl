abstract type AbstractObjective end

struct Objective{D, P, F1, F2, F3}
  domain::D
  post_processor::P
  value::F1
  gradient::F2
  hessian::F3
end

function FiniteElementContainers.create_unknowns(obj::Objective)
  return FiniteElementContainers.create_unknowns(obj.domain.static.dof)
end

function solve!(obj::Objective, solver)
  Uu = create_unknowns(obj)
  times = obj.domain.cache.times
  n = 1
  while times.current_time < times.end_time
    @info "$(repeat('=', 96))"
    @info "= Load step    $(times.current_time_step)"
    @info "= Old Time     $(times.current_time)"
    @info "= New Time     $(times.current_time + times.Δt)"
    @info "= End Time     $(times.end_time)"
    @info "= % Completete $(100.0 * (times.current_time + times.Δt) / times.end_time)"
    @info "$(repeat('=', 96))"
    n += 1
    step!(obj.domain)
    solve!(solver, obj, Uu)
    write_time(obj.post_processor, n, times.current_time)
    write_values(obj.post_processor, NodalVariable, n, "displ_x", obj.domain.cache.U[1, :])
    write_values(obj.post_processor, NodalVariable, n, "displ_y", obj.domain.cache.U[2, :])
  end
  close(obj.post_processor)
  return nothing
end

# a little dirty but whatever
gradient(obj::Objective, Uu) = obj.gradient(obj.domain, Uu)
hessian(obj::Objective, Uu) = obj.hessian(obj.domain, Uu)
value(obj::Objective, Uu) = obj.value(obj.domain, Uu)

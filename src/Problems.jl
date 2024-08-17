abstract type AbstractProblem{S, O, D, P} end

function create_unknowns(prob::P) where P <: AbstractProblem
  return create_unknowns(prob.domain)
end

function write_time(prob::P) where P <: AbstractProblem
  times = prob.domain.times
  write_time(prob.post_processor, times.current_time_step, times.current_time)
end

# TODO 
# modify to handle arbitrary names
function write_values(prob::P, ::Type{NodalVariable}) where P <: AbstractProblem
  n = prob.domain.times.current_time_step
  write_values(prob.post_processor, NodalVariable, n, "displ_x", prob.domain.U[1, :])
  write_values(prob.post_processor, NodalVariable, n, "displ_y", prob.domain.U[2, :])
end

struct QuasiStaticProblem{S, O, D, P} <: AbstractProblem{S, O, D, P}
  solver::S
  objective::O
  domain::D
  post_processor::P
end

function solve!(Uu, prob::QuasiStaticProblem)
  @unpack solver, objective, domain, post_processor = prob
  times = domain.times
  while times.current_time < (times.end_time - times.Δt)
    @info "$(repeat('=', 96))"
    @info "= Load step    $(times.current_time_step)"
    @info "= Old Time     $(times.current_time)"
    @info "= New Time     $(times.current_time + times.Δt)"
    @info "= End Time     $(times.end_time)"
    @info "= % Completete $(100.0 * (times.current_time + times.Δt) / times.end_time)"
    @info "$(repeat('=', 96))"
    step!(domain)
    solve!(
      Uu, solver, objective, domain
      # linear_solver_alg=linear_solver_alg
    )
    write_time(prob)
    write_values(prob, NodalVariable)
  end
  close(prob.post_processor)
  return nothing
end

function solve(prob::QuasiStaticProblem)
  Uu = create_unknowns(prob)
  # TODO temporary fix until we fix up Exodus.jl
  try
    solve!(Uu, prob)
  catch e
    @warn "Something bad happend so closing output files"
    close(prob.post_processor)
    throw(e)
  end
  return Uu
end

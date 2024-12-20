"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct QuasiStaticProblem{O, S, P, T} <: AbstractProblem{O, S, P, T}
  objective::O
  solver::S
  post_processor::P
  timer::T
end

function QuasiStaticProblem(inputs::Dict{Symbol, Any})
  timer = TimerOutput()
  @timeit timer "QuasiStaticProblem - setup" begin
    domain_inputs = inputs[:domain]
    obj_inputs = inputs[:objective]
    time_inputs = inputs[:time]
    solver_inputs = inputs[Symbol("nonlinear solver")]
    pp_inputs = inputs[:postprocessor]

    domain = eval(Symbol(domain_inputs[:type]))(domain_inputs)
    obj = eval(Symbol(obj_inputs[:type]))(obj_inputs, domain, timer)
    time = eval(Symbol(time_inputs[:type]))(time_inputs)
    p = ObjectiveParameters(obj, time)
    solver = eval(Symbol(solver_inputs[:type]))(solver_inputs, obj, p, timer)
    pp = eval(Symbol(pp_inputs[:type]))(pp_inputs, domain_inputs[Symbol("mesh file")])
    return QuasiStaticProblem(obj, solver, pp, timer), p
  end
end

"""
$(TYPEDSIGNATURES)
"""
function load_step_banner(::QuasiStaticProblem, times)
  @info "$(repeat('=', 96))"
  @info "= Load step    $(times.current_time_step[1])"
  @info "= Old Time     $(times.current_time[1])"
  @info "= New Time     $(times.current_time[1] + times.Δt[1])"
  @info "= End Time     $(times.end_time[1])"
  @info "= % Completete $(100.0 * (times.current_time[1] + times.Δt[1]) / times.end_time[1])"
  @info "$(repeat('=', 96))"
end

"""
$(TYPEDSIGNATURES)
"""
function solve!(prob::QuasiStaticProblem, Uu, p)
  @unpack objective, solver, post_processor = prob
  @timeit timer(prob) "QuasiStaticProblem - solve!" begin
    try
      # write initial conditions as step 1
      n = 1
      update_dirichlet_vals!(p, objective)
      update_neumann_vals!(p, objective)
      update_field_unknowns!(p.U, objective.domain, Uu)
      write_time(post_processor, 1, 0.0) # TODO init time might not be zero
      # TODO we may want more outputs
      write_fields(post_processor, p.U, 1)

      # while p.t.current_time < p.t.end_time
      while current_time(p.t) < end_time(p.t)
        load_step_banner(prob, p.t)
        @timeit timer(prob) "Load step - solve!" begin
          step!(p)
          update_dirichlet_vals!(p, objective)
          update_neumann_vals!(p, objective)
          solve!(solver, Uu, p)
        end

        # post processing TODO lots todo
        @timeit timer(prob) "Load step - postprocess" begin
          update_field_unknowns!(p.U, objective.domain, Uu)
          write_time(post_processor, n + 1, current_time(p.t))
          write_fields(post_processor, p.U, n + 1)
        end
        n = n + 1
      end
    catch e
      close(post_processor)
      throw(e)
    end
    close(post_processor)
  end

  # show(merge(merge(solver), objective.timer))
  show(prob.timer)
end

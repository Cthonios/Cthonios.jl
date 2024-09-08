abstract type AbstractProblem{O, S, P, T} end

struct QuasiStaticProblem{O, S, P, T} <: AbstractProblem{O, S, P, T}
  objective::O
  solver::S
  post_processor::P
  timer::T
end

timer(prob::QuasiStaticProblem) = prob.timer

function load_step_banner(::QuasiStaticProblem, times)
  @info "$(repeat('=', 96))"
  @info "= Load step    $(times.current_time_step)"
  @info "= Old Time     $(times.current_time)"
  @info "= New Time     $(times.current_time + times.Δt)"
  @info "= End Time     $(times.end_time)"
  @info "= % Completete $(100.0 * (times.current_time + times.Δt) / times.end_time)"
  @info "$(repeat('=', 96))"
end

function solve!(prob::QuasiStaticProblem, Uu, p)
  @unpack objective, solver, post_processor = prob
  @timeit timer(prob) "QuasiStaticProblem - solve!" begin
    try
      # write initial conditions as step 1
      n = 1
      update_dirichlet_vals!(p, objective)
      update_field_unknowns!(p.U, objective.domain, Uu)
      write_time(post_processor, 1, 0.0) # TODO init time might not be zero
      # TODO we may want more outputs
      write_fields(post_processor, p.U, 1)

      while p.t.current_time < p.t.end_time
        load_step_banner(prob, p.t)
        @timeit timer(prob) "Load step - solve!" begin
          if solver.use_warm_start
            # TODO currently stepping is done with in the warm start solver
            # figure out how to sort this out and make it  consistent across
            # interfaces
            solve!(solver.warm_start, solver.preconditioner, solver.objective, Uu, p)
          else
            step!(p)
            update_dirichlet_vals!(p, objective)
          end
          solve!(solver, Uu, p)
        end

        # post processing TODO lots todo
        @timeit timer(prob) "Load step - postprocess" begin
          update_field_unknowns!(p.U, objective.domain, Uu)
          write_time(post_processor, n + 1, p.t.current_time)
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

# exports
export QuasiStaticProblem

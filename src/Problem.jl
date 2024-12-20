struct Problem{I, S, P, T}
  integrator::I
  solver::S
  post_processor::P
  timer::T
end

function Problem(
  objective::AbstractObjective, integrator, solver_type, pp; use_warm_start=true
) # TODO add solver options
  p = ObjectiveParameters(objective, integrator)
  solver = solver_type(objective, p, objective.timer; use_warm_start=use_warm_start)
  return Problem(integrator, solver, pp, objective.timer)
end

function Problem(inputs::Dict{Symbol, Any})
  timer = TimerOutput()
  @timeit timer "Problem - setup" begin
    domain_inputs = inputs[:domain]
    obj_inputs = inputs[:objective]
    time_inputs = inputs[:time]
    solver_inputs = inputs[Symbol("nonlinear solver")]
    pp_inputs = inputs[:postprocessor]

    domain = eval(Symbol(domain_inputs[:type]))(domain_inputs)
    obj = eval(Symbol(obj_inputs[:type]))(obj_inputs, domain, timer)
    integrator = eval(Symbol(time_inputs[:type]))(time_inputs)
    p = ObjectiveParameters(obj, integrator)
    solver = eval(Symbol(solver_inputs[:type]))(solver_inputs, obj, p, timer)
    pp = eval(Symbol(pp_inputs[:type]))(pp_inputs, domain_inputs[Symbol("mesh file")])
    return Problem(integrator, solver, pp, timer), p
  end
end

timer(prob::Problem) = prob.timer

function create_unknowns_and_parameters(prob::Problem)
  int_unknowns = integrator_unknowns(prob.solver.objective, prob.integrator)
  p = ObjectiveParameters(prob.solver.objective, prob.integrator)
  return int_unknowns..., p
end

function solve!(prob::Problem, Uu, p)
  @unpack integrator, solver, post_processor = prob
  objective = solver.objective
  @timeit timer(prob) "Problem - solve!" begin
    try
      n = 1
      update_dirichlet_vals!(p, objective)
      update_neumann_vals!(p, objective)
      update_field_unknowns!(current_solution(p), objective.domain, Uu)

      write_output(post_processor, 1, objective, Uu, p)

      while current_time(integrator) < end_time(integrator)
        @timeit timer(prob) "Integrator step" begin
          # integration_step_header(integrator)
          # step_new!(p, objective)
          # solve!(solver, Uu, p)
          step!(integrator, solver, Uu, p)
        end

        write_output(post_processor, n + 1, objective, Uu, p)

        n = n + 1        
      end
    catch e
      close(post_processor)
      throw(e)
    end
    close(post_processor)
  end
  show(prob.timer)
  return nothing
end

export Problem

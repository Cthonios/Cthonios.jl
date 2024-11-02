"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct EigenProblem{O, S, P, T} <: AbstractProblem{O, S, P, T}
  objective::O
  solver::S
  post_processor::P
  timer::T
end

function solve!(prob::EigenProblem, Uu, p)
  @unpack objective, solver, post_processor = prob
  @timeit timer(prob) "EigenProblem - solve!" begin
    try
      update_dirichlet_vals!(p, objective)
      update_neumann_vals!(p, objective)
      update_field_unknowns!(p.U, objective.domain, Uu)
      results = solve!(solver, objective, Uu, p)
      display(results)
      lambdas = sqrt.(1. ./ results.λ)
      display(lambdas)

      U_preload = copy(p.U)
      for n in axes(lambdas, 1)
        # update_field_unknowns!(p.U, objective.domain, results.X[:, n])
        update_field_unknowns!(p.U, objective.domain, results.X[:, n])
        write_time(post_processor, n, lambdas[n])
        write_fields(post_processor, p.U .+ U_preload, n)
        # write_fields(post_processor, p.U, n)
      end 
    catch e
      close(post_processor)
      throw(e)
    end
    close(post_processor)
  end
  show(prob.timer)
end
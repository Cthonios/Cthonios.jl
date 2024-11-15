# TODO refactor below to have cached values for the 
# variables needed fro warmstarting

struct WarmStart{U1, U2}
  Uu_old::U1
  U_old::U2
end

function WarmStart(o::Objective)
  return WarmStart(create_unknowns(o.domain), create_fields(o.domain))
end

# solve! in EnzymeExt currently.
# TODO figure out a way to do it analytically
function solve!(::WarmStart, solver, objective, Uu, p)
  @info "Warm start"

  # need scratch arrays
  R = solver.assembler.residuals
  R .= zero(eltype(R))
  
  dR = make_zero(R)
  dUu = make_zero(Uu)
  dp = make_zero(p)

  # need to set time step
  dp.t.current_time[1] = p.t.current_time[1]
  dp.t.current_time_step[1] = p.t.current_time_step[1]
  dp.t.Δt[1] = p.t.Δt[1]
  Cthonios.step!(dp.t)

  # Cthonios.step!(dp)

  # @show dp.t.current_time
  Cthonios.update_dirichlet_vals!(dp, objective)
  Cthonios.update_neumann_vals!(dp, objective)
  # display(dp.Ubc |> sum)
  # display(dp.nbc |> sum)
  dp.Ubc .= p.Ubc .- dp.Ubc
  dp.nbc .= p.nbc .- dp.nbc
  autodiff(
    Forward, Cthonios.gradient_for_ad!,
    Duplicated(R, dR),
    Const(objective),
    Duplicated(Uu, dUu),
    Duplicated(p, dp)
  )

  # need to subtract 1 from teh time step int since it steps twice here
  p.t.current_time_step[1] = p.t.current_time_step[1] - 1

  R = dR[objective.domain.dof.unknown_dofs]

  K = Cthonios.hessian!(solver.assembler, objective, Uu, p)
  dx = K \ R
  
  Uu .= Uu .+ dx
  return nothing
end

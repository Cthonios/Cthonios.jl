struct WarmStart{R, Uu, p}
  dR::R
  dUu::Uu
  dp::p
  ΔUu
end

function WarmStart(o::AbstractObjective, p)
  dR = create_fields(o.domain)
  dUu = create_unknowns(o.domain)
  dp = make_zero(p)
  ΔUu = create_unknowns(o.domain)
  return WarmStart(dR, dUu, dp, ΔUu)
end

# TODO figure out a way to do it analytically
function solve!(warm_start::WarmStart, solver, objective, Uu, p)
  @timeit timer(objective) "Warm start - solve!" begin
    @info "Warm start"

    # need scratch arrays
    R = solver.assembler.residuals
    R .= zero(eltype(R))
    
    @unpack dR, dUu, dp, ΔUu = warm_start

    dR .= zero(eltype(dR))
    dUu .= zero(eltype(dUu))
    zero_parameters!(dp)

    # need to set time step
    dp.t.current_time[1] = p.t.current_time[1]
    dp.t.current_time_step[1] = p.t.current_time_step[1]
    dp.t.Δt[1] = p.t.Δt[1]
    Cthonios.step!(dp.t)

    # set bcs
    Cthonios.update_dirichlet_vals!(dp, objective)
    Cthonios.update_neumann_vals!(dp, objective)
    dp.Ubc .= p.Ubc .- dp.Ubc
    dp.nbc .= p.nbc .- dp.nbc

    # TODO try to make this use finite difference
    # instead of having to rely on Enzyme since
    # Enyzme takes forever to compile/precompile
    @timeit timer(objective) "WarmStart - AD" begin
      autodiff(
        Forward, Cthonios.gradient_for_ad!,
        Duplicated(R, dR),
        Const(objective),
        Duplicated(Uu, dUu),
        Duplicated(p, dp)
      )
    end

    # need to subtract 1 from teh time step int since it steps twice here
    p.t.current_time_step[1] = p.t.current_time_step[1] - 1

    R = dR[objective.domain.dof.unknown_dofs]

    K = Cthonios.hessian!(solver.assembler, objective, Uu, p)
    @timeit timer(objective) "WarmStart - solve" begin
      # ΔUu .= zero(eltype(ΔUu))
      # IterativeSolvers.gmres!(ΔUu, K, R; verbose=true)
      # ldiv!(ΔUu, K, R)
      ΔUu .= K \ R
    end
    
    Uu .= Uu .+ ΔUu
    return nothing
  end
end

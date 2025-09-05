struct WarmStart{RT, Uu, p}#, S}
  R::RT
  dR::RT
  dUu::Uu
  dp::p
  ΔUu::Uu
  # solver::S
end

function WarmStart(o::AbstractObjectiveCache, p)
  R = create_field(o)
  dR = create_field(o)
  dUu = create_unknowns(o)
  dp = make_zero(p)
  ΔUu = create_unknowns(o)
  # solver = GmresSolver(length(dUu), length(dUu), length(dUu), typeof(dUu))
  return WarmStart(R, dR, dUu, dp, ΔUu)#, solver)
end

# TODO figure out a way to do it analytically
function solve!(warm_start::WarmStart, objective, Uu, p; verbose=false)
  @timeit objective.timer "Warm start - solve!" begin
    if verbose
      @info "Warm start"
    end

    assembler = objective.sim_cache.assembler

    (; R, dR, dUu, dp, ΔUu) = warm_start
    fill!(R, zero(eltype(R)))
    fill!(dR, zero(eltype(dR)))
    fill!(dUu, zero(eltype(dUu)))

    dp = make_zero(dp)
    # remake_zero!(dp)

    # need to set time step
    fill!(dp.times.time_current, sum(p.times.time_current))
    fill!(dp.times.Δt, sum(p.times.Δt))
    FiniteElementContainers.update_time!(dp)

    # set bcs
    FiniteElementContainers.update_bc_values!(
      dp.dirichlet_bcs, dp.dirichlet_bc_funcs, p.h1_coords, sum(dp.times.time_current)
    )
    FiniteElementContainers.update_bc_values!(
      dp.neumann_bcs, dp.neumann_bc_funcs, p.h1_coords, sum(dp.times.time_current)
    )

    # TODO needs to be updated for GPUs
    for (bc, dbc) in zip(values(p.dirichlet_bcs), values(dp.dirichlet_bcs))
      dbc.vals .= bc.vals .- dbc.vals
    end

    for (bc, dbc) in zip(values(p.neumann_bcs), values(dp.neumann_bcs))
      dbc.vals .= bc.vals .- dbc.vals
    end

    @timeit objective.timer "WarmStart - AD" begin
      autodiff(
        Forward, assemble_vector_for_ad!,
        Duplicated(R, dR),
        Const(assembler),
        Duplicated(Uu, dUu),
        Duplicated(p, dp),
        Const(residual)
      )
    end

    # need to subtract 1 from teh time step int since it steps twice here
    # p.t.current_time_step[1] = p.t.current_time_step[1] - 1

    # TODO needs to be updated for GPUs
    # R = dR[objective.sim_cache.assembler.dof.H1_unknown_dofs]
    R = dR[objective.sim_cache.assembler.dof.unknown_dofs]
    # R = residual(assembler)

    # K = Cthonios.hessian!(solver.assembler, objective, Uu, p)
    assemble_stiffness!(assembler, objective.objective.hessian_u, Uu, p)
    # K = hessian(objective, Uu, p)
    # TODO below will fail for dynamics
    K = stiffness(assembler)
    @timeit objective.timer "WarmStart - solve" begin
      ΔUu .= K \ R
      # ldiv!(ΔUu, K, R)
      # ldiv!(K, R)
      # copyto!(ΔUu, R)
      # ldiv!(ΔUu, K, R)
      # Krylov.solve!(warm_start.solver, K, R)
      # copyto!(ΔUu, Krylov.solution(warm_start.solver))
    end
    
    # Uu .= Uu .+ ΔUu
    Uu .+= ΔUu
    return nothing
  end
end

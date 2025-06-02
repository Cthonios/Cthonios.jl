struct WarmStart{R, Uu, p}
  dR::R
  dUu::Uu
  dp::p
  ΔUu
end

function WarmStart(o::AbstractObjective, p)
  dR = create_field(o)
  dUu = create_unknowns(o)
  dp = make_zero(p)
  ΔUu = create_unknowns(o)
  return WarmStart(dR, dUu, dp, ΔUu)
end

function assemble_vector_for_ad!(storage, assembler, Uu, p)
  fill!(storage, zero(eltype(storage)))
  fspace = FiniteElementContainers.function_space(assembler, H1Field)
  t = FiniteElementContainers.current_time(p.times)
  Δt = FiniteElementContainers.time_step(p.times)
  update_bcs!(p)
  update_field_unknowns!(p.h1_field, assembler.dof, Uu)
  for (b, (conns, block_physics, state_old, state_new, props)) in enumerate(zip(
    values(fspace.elem_conns), 
    values(p.physics),
    values(p.state_old), values(p.state_new),
    values(p.properties)
  ))
    ref_fe = values(fspace.ref_fes)[b]
    backend = FiniteElementContainers._check_backends(assembler, p.h1_field, p.h1_coords, state_old, state_new, conns)
    FiniteElementContainers._assemble_block_vector!(
      storage, block_physics, ref_fe, 
      p.h1_field, p.h1_coords, state_old, state_new, props, t, Δt,
      conns, b, residual,
      backend
    )
  end
end

# TODO figure out a way to do it analytically
function solve!(warm_start::WarmStart, objective, Uu, p)
  @timeit objective.timer "Warm start - solve!" begin
    @info "Warm start"
    assembler = objective.assembler
    # need scratch arrays
    # R = assembler.residual_unknowns
    R = create_field(assembler, H1Field)
    R .= zero(eltype(R))
    
    # @unpack dR, dUu, dp, ΔUu = warm_start
    dR, dUu, dp_temp, ΔUu = warm_start.dR, warm_start.dUu, warm_start.dp, warm_start.ΔUu

    dR .= zero(eltype(dR))
    dUu .= zero(eltype(dUu))
    # zero_parameters!(dp)
    # TODO major inefficiency here
    dp = make_zero(dp_temp)

    # need to set time step
    dp.times.time_current[1] = p.times.time_current[1]
    # dp.times.Δt[1] = p.time.current_time_step[1]
    dp.times.Δt[1] = p.times.Δt[1]
    # Cthonios.step!(dp.t)
    FiniteElementContainers.update_time!(dp)

    # set bcs
    # Cthonios.update_dirichlet_vals!(dp, objective)
    # Cthonios.update_neumann_vals!(dp, objective)
    # dp.Ubc .= p.Ubc .- dp.Ubc
    # dp.nbc .= p.nbc .- dp.nbc
    FiniteElementContainers.update_bc_values!(
      dp.dirichlet_bcs, dp.dirichlet_bc_funcs, p.h1_coords, dp.times.time_current[1]
    )

    for (bc, dbc) in zip(values(p.dirichlet_bcs), values(dp.dirichlet_bcs))
      dbc.vals .= bc.vals .- dbc.vals
    end

    # @assert false
    # TODO try to make this use finite difference
    # instead of having to rely on Enzyme since
    # Enyzme takes forever to compile/precompile
    @timeit objective.timer "WarmStart - AD" begin
      autodiff(
        # Forward, Cthonios.gradient_for_ad!,
        # Duplicated(R, dR),
        # Const(objective),
        # Forward, Cthonios.gradient,
        # Duplicated(objective, dobjective),
        # Const(objective),
        Forward, assemble_vector_for_ad!,
        Duplicated(R, dR),
        Const(assembler),
        Duplicated(Uu, dUu),
        Duplicated(p, dp)
      )
    end

    # need to subtract 1 from teh time step int since it steps twice here
    # p.t.current_time_step[1] = p.t.current_time_step[1] - 1

    R = dR[objective.assembler.dof.H1_unknown_dofs]
    # R = residual(assembler)

    # K = Cthonios.hessian!(solver.assembler, objective, Uu, p)
    K = hessian(objective, Uu, p)
    @timeit objective.timer "WarmStart - solve" begin
      # ΔUu .= zero(eltype(ΔUu))
      # IterativeSolvers.gmres!(ΔUu, K, R; verbose=true)
      # ldiv!(ΔUu, K, R)
      ΔUu .= K \ R
    end
    
    Uu .= Uu .+ ΔUu
    return nothing
  end
end

struct WarmStart{RT, Uu, p}#, S}
  R::RT
  dR::RT
  dUu::Uu
  dp::p
  ΔUu::Uu
  # solver::S
  timer::TimerOutput
end

function WarmStart(o, p, timer)
  @timeit timer "WarmStart - setup" begin
    # Uu = o.solution
    # p = o.parameters
    R = make_zero(o.gradient)
    dR = make_zero(o.gradient)
    # dUu = make_zero(o.solution_old.data)
    dUu = create_unknowns(o)
    dp = make_zero(p)
    # ΔUu = make_zero(o.solution_old.data)
    ΔUu = create_unknowns(o)
    # solver = GmresSolver(length(dUu), length(dUu), length(dUu), typeof(dUu))

    return WarmStart(R, dR, dUu, dp, ΔUu, timer)#, solver)
  end
end

# TODO figure out a way to do it analytically
function solve!(
  warm_start::WarmStart, objective, Uu, p;
  P=I,
  verbose=false
)
  @timeit warm_start.timer "Warm start - solve!" begin
    if verbose
      @info "Warm start"
    end

    asm = assembler(objective)

    (; R, dR, dUu, dp, ΔUu) = warm_start
    fill!(R, zero(eltype(R)))
    fill!(dR, zero(eltype(dR)))
    fill!(dUu, zero(eltype(dUu)))

    # set params to zeros
    # remake_zero!(dp)
    dp = make_zero(dp)

    # set bcs
    @timeit warm_start.timer "warm start - bc updates" begin
      # need to set time step
      fill!(dp.times.time_current, sum(p.times.time_current))
      fill!(dp.times.Δt, sum(p.times.Δt))
      FiniteElementContainers.update_time!(dp)

      # update values
      FiniteElementContainers.update_bc_values!(
        dp.dirichlet_bcs, p.h1_coords, sum(dp.times.time_current)
      )
      FiniteElementContainers.update_bc_values!(
        dp.neumann_bcs, p.h1_coords, sum(dp.times.time_current)
      )

      # TODO needs to be updated for GPUs
      for (bc, dbc) in zip(values(p.dirichlet_bcs.bc_caches), values(dp.dirichlet_bcs.bc_caches))
        dbc.vals .= bc.vals .- dbc.vals
      end

      for (bc, dbc) in zip(values(p.neumann_bcs.bc_caches), values(dp.neumann_bcs.bc_caches))
        dbc.vals .= bc.vals .- dbc.vals
      end
    end

    @timeit warm_start.timer "WarmStart - AD" begin
      # autodiff(
      #   Forward, assemble_vector_for_ad!,
      #   Duplicated(R, dR),
      #   Const(asm),
      #   Duplicated(Uu, dUu),
      #   Duplicated(p, dp),
      #   Const(residual)
      # )
      autodiff(
        Forward, 
        assemble_vector!,
        Duplicated(R, dR),
        Const(asm.dof),
        Const(residual),
        Duplicated(Uu, dUu),
        Duplicated(p, dp)
      )

      # TODO needs to be updated for GPUs
      if FiniteElementContainers._is_condensed(asm.dof)
        # do nothing
        R = dR.data
      else
        R = dR[asm.dof.unknown_dofs]
      end
    end

    # TODO below will fail for dynamics
    @timeit warm_start.timer "WarmStart - assembly" begin
      assemble_stiffness!(asm, objective.objective.hessian_u, Uu, p)    
      K = stiffness(asm)
    end

    # @timeit warm_start.timer "WarmStart - factorize" begin
    #   P = warm_start.P
    #   lu!(P, K)
    # end

    @timeit warm_start.timer "WarmStart - solve" begin
      ΔUu .= K \ R

      # Krylov.solve!(CgSolver(), K, R)
      # P = lu(K)
      # ΔUu, stats = Krylov.cg(K |> Symmetric, R; verbose=1, ldiv=true)
      # K = stiffness(asm)
      # ΔUu, state = Krylov.gmres(K, R; M=P, ldiv=true)
    end
    
    # Uu .= Uu .+ ΔUu
    Uu .+= ΔUu
    return nothing
  end
end

"""
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct TrustRegionSolverSettings #<: NonlinearSolverSettings
  t_1::Float64                                  = 0.25
  t_2::Float64                                  = 1.75
  η_1::Float64                                  = 1e-10
  η_2::Float64                                  = 0.1
  η_3::Float64                                  = 0.5
  max_trust_iters::Int64                        = 100
  tol::Float64                                  = 1.e-8
  max_cg_iters::Int64                           = 50
  max_cumulative_cg_iters::Int64                = 1000
  cg_tol::Float64                               = 0.2e-9
  cg_inexact_solve_ratio::Float64               = 1e-5
  tr_size::Float64                              = 2.0
  min_tr_size::Float64                           = 1e-8
  check_stability::Bool                         = false
  use_preconditioned_inner_product_for_cg::Bool = true
  use_incremental_objective::Bool               = false
  debug_info::Bool                              = true
  over_iter::Int64                              = 0
end

function Base.show(io::IO, settings::TrustRegionSolverSettings)
  print(io, "  TrustRegionSolverSettings\n")
  for name in fieldnames(typeof(settings))
    print(io, "        ", rpad(name, 39), " = ", getfield(settings, name), "\n")
  end
  print(io, "\n")
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct TrustRegionSolver{
  L,
  O,
  U <: AbstractVector,
  W,
  T <: TimerOutput,
  F <: H1Field,
  S <: TrustRegionSolverSettings
} #<: AbstractNonlinearSolver{L, O, U, W, T}
  preconditioner::L
  objective_cache::O
  ΔUu::U
  warm_start::W
  timer::T
  o::U
  g::F
  Hv::F
  cauchy_point::U
  q_newton_point::U
  d::U
  y_scratch_1::U
  y_scratch_2::U
  y_scratch_3::U
  y_scratch_4::U
  use_warm_start::Bool
  verbose::Bool
  settings::S
end

"""
$(TYPEDSIGNATURES)
TODO figure out which scratch arrays can be nixed
"""
function TrustRegionSolver(
  objective; 
  preconditioner=CholeskyPreconditioner,
  use_warm_start=true,
  verbose=true,
  settings=TrustRegionSolverSettings()
)
  timer = TimerOutput()

  @timeit timer "TrustRegionSolver - setup" begin
    # domain = objective.domain
    # settings       = TrustRegionSolverSettings() # TODO add non-defaults
    # TODO eventually write a custom linear solver for this one
    p              = objective.parameters
    precond        = preconditioner(objective, timer)
    ΔUu            = create_unknowns(objective)
    # TODO
    warm_start     = WarmStart(objective)
    # warm_start     = nothing
    o              = zeros(1)
    g              = create_field(objective) # gradient
    Hv             = create_field(objective) # hessian-vector product
    cauchy_point   = create_unknowns(objective)
    q_newton_point = create_unknowns(objective)
    d              = create_unknowns(objective)
    y_scratch_1    = create_unknowns(objective)
    y_scratch_2    = create_unknowns(objective)
    y_scratch_3    = create_unknowns(objective)
    y_scratch_4    = create_unknowns(objective)
  end
  return TrustRegionSolver(
    precond, objective,
    ΔUu,
    warm_start,
    timer,
    o, g, Hv,
    cauchy_point, q_newton_point, d,
    y_scratch_1, y_scratch_2, y_scratch_3, y_scratch_4,
    use_warm_start,
    verbose,
    settings
  )
end

function TrustRegionSolver(
  inputs::Dict{Symbol, Any},
  objective::AbstractObjectiveCache
)
  preconditioner = eval(Symbol(inputs[:preconditioner][:type]))
  warm_start = inputs[Symbol("warm start")]
  settings = TrustRegionSolverSettings()

  for key in fieldnames(typeof(settings))
    if haskey(inputs, key)
      setfield!(settings, key, inputs[key])
    end
  end

  return TrustRegionSolver(objective; preconditioner=preconditioner, use_warm_start=warm_start, settings=settings)
end

timer(solver::TrustRegionSolver) = solver.timer

negCurveString = "neg curve"
boundaryString = "boundary"
interiorString = "interior"

function is_converged(solver::TrustRegionSolver, objective, x, realO, modelO, realRes, modelRes, cgIters, trSize, settings)
  gg = dot(realRes, realRes)
  if gg < settings.tol^2
    modelResNorm = norm(modelRes)
    realResNorm = sqrt(gg)
    print_min_banner(
      solver,
      realO, modelO,
      realResNorm,
      modelResNorm,
      cgIters,
      trSize,
      interiorString,
      true,
      settings
    )

    if solver.verbose
      @info "Converged"
      println("") # a bit of output formatting
    end

    if settings.check_stability
      @assert false "hook this up"
      # objective.check_stability(x)
    end
    
        
    return true
  end
  return false
end

function is_on_boundary(stepType)
  return stepType==boundaryString || stepType==negCurveString
end

function cg_inner_products_preconditioned(alpha, beta, zd, dd, rPr, z, d)
  # recurrence formulas from Gould et al. doi:10.1137/S1052623497322735
  zd = beta * (zd + alpha * dd)
  dd = rPr + beta * beta * dd
  return zd, dd
end

function cg_inner_products_unpreconditioned(alpha, beta, zd, dd, rPr, z, d)
  zd = dot(z, d)
  dd = dot(d, d)
  return zd, dd
end

function preconditioned_project_to_boundary(z, d, trSize, zz, mult_by_approx_hessian)
  # find tau s.t. (z + tau*d)^2 = trSize^2
  Pd = mult_by_approx_hessian(d)
  dd = dot(d, Pd)
  zd = dot(z, Pd)
  tau = (sqrt((trSize^2 - zz) * dd + zd^2) - zd) / dd
  # is it ever better to choose the - sqrt() branch?
  return z + tau * d
end

function print_banner(solver::TrustRegionSolver)
  if solver.verbose
    @info "=============================================================================================================================="
    @info "Objective        Model Objective  Residual        Model Residual      CG    TR size"
    @info "=============================================================================================================================="
  end
end

function print_min_banner(
  solver::TrustRegionSolver,
  objective, modelObjective, res, modelRes, 
  cgIters, trSize, onBoundary, willAccept, settings
)
  if settings.debug_info == false
    return
  end
  
  if willAccept
    will_accept = true
  else
    will_accept = false
  end

  if solver.verbose
    str = @sprintf "%1.6e    %1.6e    %1.6e    %1.6e    %6i    %1.6e    %s    %s" objective modelObjective res modelRes cgIters trSize onBoundary will_accept
    @info str
  end
end

function project_to_boundary_with_coefs(z, d, trSize, zz, zd, dd)
  # find tau s.t. (z + tau*d)^2 = trSize^2
  tau = (sqrt((trSize^2 - zz) * dd + zd^2) - zd) / dd
  return z + tau * d
end

function update_step_length_squared(alpha, zz, zd, dd)
  return zz + 2 * alpha * zd + alpha * alpha * dd
end

function solve_trust_region_minimization(solver, x, r, hess_vec_func, P, trSize, settings)
  # minimize r@z + 0.5*z@J@z
  z = 0. * x
  zz = 0.

  cgInexactRelTol = settings.cg_inexact_solve_ratio
  cgTolSquared = max(settings.cg_tol^2, cgInexactRelTol*cgInexactRelTol * dot(r, r))
  if dot(r, r) < cgTolSquared
    return z, z, interiorString, 0
  end
  
  Pr = P \ r
  # Pr = solver.y_scratch_4
  # ldiv!(Pr, P, r)
  d = -Pr
  cauchyP = d
  rPr = dot(r, Pr)

  zz = 0.0
  zd = 0.0
  if settings.use_preconditioned_inner_product_for_cg
    dd = rPr
    cg_inner_products = cg_inner_products_preconditioned
  else
    dd = dot(d, d)
    cg_inner_products = cg_inner_products_unpreconditioned
  end

  for i in 1:settings.max_cg_iters
    curvature = dot(d, hess_vec_func(d))
    alpha = rPr / curvature
        
    zNp1 = z + alpha*d
    zzNp1 = update_step_length_squared(alpha, zz, zd, dd)

    if curvature <= 0
      zOut = project_to_boundary_with_coefs(z, d, trSize,
                                            zz, zd, dd)
      return zOut, cauchyP, negCurveString, i+1
    end
    if zzNp1 > trSize^2
      zOut = project_to_boundary_with_coefs(z, d, trSize,
                                            zz, zd, dd)
      return zOut, cauchyP, boundaryString, i+1
    end
    z = zNp1 # z + alpha*d
        
    r += alpha * hess_vec_func(d)
    Pr = P \ r
    rPrNp1 = dot(r, Pr)
      
    if dot(r, r) < cgTolSquared
      return z, cauchyP, interiorString, i+1
    end
    beta = rPrNp1 / rPr
    rPr = rPrNp1
    d = -Pr + beta*d

    zz = zzNp1
    zd, dd = cg_inner_products(alpha, beta, zd, dd, rPr, z, d)
  end 
        
  return z, cauchyP, interiorString * "_", settings.max_cg_iters
end

function dogleg_step(cp, newtonP, trSize, mat_mul)
  cc = dot(cp, mat_mul(cp))
  nn = dot(newtonP, mat_mul(newtonP))
  tt = trSize*trSize
    
  if cc >= tt #return cauchy point if it extends outside the tr
    #print('cp on boundary')
    return cp * sqrt(tt/cc)
  end

  if cc > nn # return cauchy point?  seems the preconditioner was not accurate?
    @info "cp outside newton, preconditioner likely inaccurate"
    return cp
  end

  if nn > tt # on the dogleg (we have nn >= cc, and tt >= cc)
      #print('dogleg')
    return preconditioned_project_to_boundary(
      cp,
      newtonP-cp,
      trSize,
      cc,
      mat_mul
    )
  end
  #print('quasi-newton step')
  return newtonP

end

function solve!(solver::TrustRegionSolver, Uu, p)
  @timeit timer(solver) "TrustRegionSolver - solve!" begin
    if solver.use_warm_start
      @timeit timer(solver) "TrustRegionSolver - warm start" begin
        # solve!(solver.warm_start, solver.preconditioner, solver.objective_cache, Uu, p)
        solve!(solver.warm_start, solver.objective_cache, Uu, p; verbose=solver.verbose)
      end
    end

    # unpack some solver settings
    settings = solver.settings
    tr_size = solver.settings.tr_size

    # set up some stuff
    # TODO clean up later
    x = Uu

    # saving initial objective and gradient
    o = value(solver.objective_cache, x, p)
    g = gradient(solver.objective_cache, x, p)
    g_norm = norm(g)

    if solver.verbose
      @info @sprintf "Initial objective = %1.6e" o[1]
      @info @sprintf "Initial residual  = %1.6e" g_norm
    end

    # preconditioner
    update_preconditioner!(solver.preconditioner, solver.objective_cache, x, p)
    P = solver.preconditioner.preconditioner
    # this could potentially return an unstable solution
    if is_converged(solver, value, x, 0.0, 0.0, g, g, 0, tr_size, settings)
      # if callback
      #   # callback(x, objective)
      #   @assert false
      # end
      # return x, True
      if solver.verbose
        @show "converged"
      end
      return nothing
    end

    cumulativeCgIters=0
    print_banner(solver)
    for i in 1:settings.max_trust_iters
      if i > 1 && i % 25 == 0
        print_banner(solver)
      end

      # setup some local closures to this iteration
      if settings.use_incremental_objective
        @assert false "not implemented yet"
      else
        increment_objective = d -> value(solver.objective_cache, x + d, p) - o
      end

      hess_vec_func = v -> hvp(solver.objective_cache, x, p, v)
      # TODO need to fix below
      # K = hessian!(solver.preconditioner.assembler, solver.objective_cache, x, p)
      # mult_by_approx_hessian = v -> (K + 0.0 * I) * v
      # mult_by_approx_hessian = v -> K \ v
      # mult_by_approx_hessian = v -> hvp(solver.objective_cache, x, p, v)
      mult_by_approx_hessian = v -> v

      # calculate cauchy point
      @timeit timer(solver) "TrustRegionSolver - cauchy point" begin
        gKg = dot(g, hess_vec_func(g))
        if gKg > 0
          alpha = -dot(g, g) / gKg
          cauchyPoint = alpha * g
          cauchyPointNormSquared = dot(cauchyPoint, mult_by_approx_hessian(cauchyPoint))
        else
          cauchyPoint =  -g * (tr_size / sqrt.(dot(g, mult_by_approx_hessian(g))))
          cauchyPointNormSquared = tr_size * tr_size
          @info "negative curvature unpreconditioned cauchy point direction found."
        end
      end

      # solve minimize
      if cauchyPointNormSquared >= tr_size*tr_size
        @info "unpreconditioned gradient cauchy point outside trust region at dist = $(sqrt.(cauchyPointNormSquared))"
        cauchyPoint *= (tr_size / sqrt(cauchyPointNormSquared))
        cauchyPointNormSquared = tr_size * tr_size
        qNewtonPoint = cauchyPoint
        stepType = boundaryString
        cgIters = 1
      else
        @timeit timer(solver) "TrustRegionSolver - minimization" begin
          qNewtonPoint, _, stepType, cgIters = 
              solve_trust_region_minimization(
                solver,
                x, g,
                hess_vec_func,
                # objective.apply_precond,
                P,
                tr_size, settings
              )
        end
      end

      # line 411 in optimism EquationSolver next below
      cumulativeCgIters += cgIters
      trSizeUsed = tr_size
      happyAboutTrSize = false

      while !happyAboutTrSize
        d = dogleg_step(cauchyPoint, qNewtonPoint, tr_size, mult_by_approx_hessian)
        Jd = hess_vec_func(d)
        dJd = dot(d, Jd)
        modelObjective = dot(g, d) + 0.5 * dJd
        
        y = x + d
        realObjective = increment_objective(d)
        gy = gradient(solver.objective_cache, y, p)
        
        if is_converged(
          solver,
          solver.objective_cache, y, realObjective, modelObjective,
          gy, g + Jd, cgIters, trSizeUsed, settings
        )
          # if callback
          #   @assert false
          #   callback(y, objective)
          # end
          # return y, true
          Uu .= y
          return nothing
        end

        modelImprove = -modelObjective
        realImprove = -realObjective

        rho = realImprove / modelImprove

        if modelObjective > 0
          @info "Found a positive model objective increase.  Debug if you see this."
          @info "modelObjective = $modelObjective"
          @info "realObjective = $realObjective"
          rho = realImprove / -modelImprove
        end

        if !(rho >= settings.η_2)  # write it this way to handle NaNs
          tr_size *= settings.t_1
        elseif rho > settings.η_3 && is_on_boundary(stepType)
          tr_size *= settings.t_2
        end

        modelRes = g + Jd
        modelResNorm = norm(modelRes)
        realResNorm = norm(gy)

        willAccept = rho >= settings.η_1 || 
                     (rho >= -0 && realResNorm <= g_norm)
        print_min_banner(
          solver,
          realObjective, modelObjective,
          realResNorm, modelResNorm,
          cgIters, trSizeUsed, stepType,
          willAccept,
          settings
        )

        if willAccept
          x = y
          g = gy
          o = value(solver.objective_cache, x, p)
          g_norm = realResNorm
          triedNewPrecond = false
          happyAboutTrSize = true
          # if callback

          # end
        else
          # set these for output
          # trust region will continue to strink until we find a solution on the boundary
          stepType=boundaryString
          cgIters = 0
        end

        if cgIters >= settings.max_cg_iters || cumulativeCgIters >= settings.max_cumulative_cg_iters
          # objective.update_precond(x)
          update_preconditioner!(solver.preconditioner, solver.objective_cache, x, p; verbose=solver.verbose)
          P = solver.preconditioner.preconditioner
          cumulativeCgIters=0
        end

        trSizeUsed = tr_size

        if tr_size < settings.min_tr_size

          if !triedNewPrecond
            @info "The trust region is too small, updating precond and trying again."
            update_preconditioner!(solver.preconditioner, solver.objective_cache, x, p)
            P = solver.preconditioner.preconditioner
            cumulativeCgIters = 0
            triedNewPrecond = true
            happyAboutTrSize = true
            tr_size = settings.tr_size                    
          else
            @info "The trust region is still too small.  Accepting, but be careful."
            # if callback: callback(x, objective)
            # return x, False
            return nothing
          end
        end
        # @assert false
      end
    end
  end

  @assert false "reached maximum tr iterations"
end
Base.@kwdef mutable struct TrustRegionSolverSettings <: NonlinearSolverSettings
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
    print(io, "    ", rpad(name, 39), " = ", getfield(settings, name), "\n")
  end
  print(io, "\n")
end

struct TrustRegionSolver{
  S <: TrustRegionSolverSettings,
  L <: LinearSolver,
  V <: AbstractVector
} <: NonlinearSolver{S, L, V}
  settings::S
  linear_solver::L
  Uu::V
  ΔUu::V
end

function TrustRegionSolver(input_settings::D, domain::QuasiStaticDomain) where D <: Dict
  settings      = TrustRegionSolverSettings() # TODO add non-defaults
  linear_solver = LinearSolver(input_settings["linear solver"], domain)
  Uu            = create_unknowns(domain)
  ΔUu           = create_unknowns(domain)
  return TrustRegionSolver(settings, linear_solver, Uu, ΔUu)
end

function Base.show(io::IO, solver::TrustRegionSolver)
  print(io, "TrustRegionSolver\n", 
        "  Settings      = $(solver.settings)\n",
        "  Linear solver = $(solver.linear_solver)\n")
end

function logger(
  ::TrustRegionSolver,
  real_obj, model_obj,
  real_res, model_res,
  cg_iters, tr_size, accept
)

  str =  @sprintf "O = %1.6e  O_M = %1.6e  ||R|| = %1.6e  ||R_R|| = %1.6e  CG = %3i  TR = %1.6e accept = " real_obj model_obj real_res model_res cg_iters tr_size
  str *= "$accept"
  @info str
end

function is_converged(
  solver::TrustRegionSolver,
  x, 
  real_obj, model_obj,
  real_res, model_res,
  cg_iters, tr_size
)

  gg = dot(real_res, real_res)
  # @show gg
  # @show solver.settings.tol^2
  if gg < solver.settings.tol^2
    model_res_norm = norm(model_res)
    real_res_norm  = sqrt(gg)
    # banner(real_obj, model_obj, real_res_norm, model_res_norm, cg_iters, tr_size, :interior)
    logger(solver, real_obj, model_obj, real_res_norm, model_res_norm, cg_iters, tr_size, true)
    # TODO need logger here

    if solver.settings.check_stability
      @assert false "not supported yet"
    end

    return true
  end

  return false
end

"""
minimize r * z + 0.5 * z * J * z
"""
function minimize_trust_region_sub_problem(
  solver::TrustRegionSolver,
  x::V, r::V, K, 
  tr_size::Float64
) where V <: AbstractVector

  z = zeros(eltype(x), length(x))

  cg_tol_squared = max(
    solver.settings.cg_inexact_solve_ratio^2 * dot(r, r),
    solver.settings.cg_tol^2
  )

  # unpack preonditioner # TODO I think it needs to be cholesky right now
  P = solver.linear_solver.solver_cache.Pl
  Pr = P \ r
  # Pr = P * r
  d = -Pr
  cauchy_point = zeros(eltype(d), length(d))
  rPr = dot(r, Pr)

  zz = zero(eltype(x))
  zd = zero(eltype(x))

  dd = rPr

  for i in 1:solver.settings.max_cg_iters
    # @info "cg iter $i"

    curvature = dot(d, K * d)
    α = rPr / curvature

    zNp1 = z + α * d
    zzNp1 = zz + 2.0 * α * zd + α * α * dd

    # TODO
    if curvature <= zero(eltype(x))
      # @assert false "Edge case not handled yet - negative curvature"
      τ = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
      z_out = z + τ * d
      return z_out, cauchy_point, :negative_curvature, i
    end

    # TODO
    if zzNp1 > tr_size^2
      # @assert false "other edge case not handled"
      # TODO projec to boundary
      τ = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
      z_out = z + τ * d
      return z_out, cauchy_point, :boundary, i
    end

    z = zNp1

    r = r + α * K * d
    Pr = P \ r
    # Pr = P * r
    rPrNp1 = dot(r, Pr)
    
    if dot(r, r) < cg_tol_squared
      return z, cauchy_point, :interior, i
    end

    β = rPrNp1 / rPr
    rPr = rPrNp1
    d = -Pr + β * d

    zz = zzNp1
    zd = β * (zd + α * dd)
    dd = rPr * β * β * dd
  end
  return z, cauchy_point, :interior, solver.settings.max_cg_iters
end

function dog_leg_step(cauchy_point, q_newton_point, tr_size, P)
  # cc = dot(cauchy_point, K * cauchy_point)
  # nn = dot(q_newton_point, K * q_newton_point)
  cc = dot(cauchy_point, P \ cauchy_point)
  nn = dot(q_newton_point, P \ q_newton_point)
  tt = tr_size * tr_size

  # return cauchy point if it extends outside the tr
  if cc >= tt
    return cauchy_point * sqrt(tt / cc)
  end

  # return cauchy point? seems the preconditioner was not accurate
  if cc > nn
    @warn "cp outside newton, preconditioner likely inaccurate"
    return cauchy_point
  end

  # on the dogleg
  if nn > tt
    # return preconditioner_project_to_boundary(cauchy_point, q_newton_point - cauchy_point, tr_size, cc, K)
    # project to boundary
    # Pd = K * (q_newton_point - cauchy_point)
    Pd = P \ (q_newton_point - cauchy_point)
    dd = dot(q_newton_point - cauchy_point, Pd)
    zd = dot(cauchy_point, Pd)
    τ  = (sqrt((tr_size^2 - cc) * dd + zd^2) - zd) / dd
    return cauchy_point + τ * (q_newton_point - cauchy_point)
  end

  return q_newton_point
end

function solve!(solver::TrustRegionSolver, domain::QuasiStaticDomain)
  # unpack cached arrays from solver
  Uu = solver.Uu
  K = stiffness(domain, Uu)
  P = solver.linear_solver.solver_cache.Pl
  # unpack some solver settings
  tr_size = solver.settings.tr_size

  # calculate initial objective and residual
  o = energy(domain, Uu)
  g = residual(domain, Uu)
  o_init = o
  g_norm_init = norm(g)
  g_norm = g_norm_init
  @info @sprintf "Initial objective = %1.6e" o_init
  @info @sprintf "Initial residual  = %1.6e" g_norm_init

  # begin loop over tr iters
  cumulative_cg_iters = 0
  for n in 1:solver.settings.max_trust_iters
    # check for convergence at beginning of iteration
    if is_converged(solver, Uu, 0.0, 0.0, g, g, 0, tr_size)
      solver.ΔUu .= solver.Uu - Uu
      solver.Uu .= Uu
    end

    objective = d -> energy(domain, Uu + d) - o
    K = stiffness(domain, Uu)

    # check for negative curvature
    gKg = dot(g, K * g)
    if gKg > 0
      α = -dot(g, g) / gKg
      cauchy_point = α * g
      cauchy_point_norm_squared = dot(cauchy_point, K * cauchy_point) # should eventually multiply by approx hessian TODO
    else
      cauchy_point = -g * (tr_size / sqrt(dot(g, K * g))) # TODO another place to mult_by_approx_hessian
      cauchy_point_norm_squared = tr_size * tr_size
      @info "negative curavture unpreconditioned cauchy point direction found."
    end

    # check if outside trust region
    if cauchy_point_norm_squared >= tr_size * tr_size
      @info "unpreconditioned gradient cauchy point outside trust region at dist = $(sqrt(cauchy_point_norm_squared))."
      cauchy_point = (tr_size / sqrt(cauchy_point_norm_squared)) * cauchy_point
      cauchy_point_norm_squared = tr_size * tr_size
      q_newton_point = cauchy_point
      step_type = :boundary
      n_cg_iters = 1
    else
      q_newton_point, _, step_type, n_cg_iters = 
      minimize_trust_region_sub_problem(solver, Uu, g, K, tr_size)
      # q_newton_point = IterativeSolvers.cg(K, -g)
      # n_cg_iters = 1
      step_type = :cg
    end

    cumulative_cg_iters += n_cg_iters

    tr_size_used = tr_size
    happy = false
    while !happy
      # take a dogleg step
      d = dog_leg_step(cauchy_point, q_newton_point, tr_size, P)

      Jd = K * d
      dJd = dot(d, Jd)
      model_objective = dot(g, d) + 0.5 * dJd

      y = Uu + d
      real_objective = objective(d)
        
      gy = residual(domain, y)

      if is_converged(solver, y, real_objective, model_objective, gy, g + Jd, n_cg_iters, tr_size_used)
        @show "converged!"
        solver.ΔUu .= solver.Uu - y
        solver.Uu .= y
        return
      end

      model_improve = -model_objective
      real_improve = -real_objective
      ρ = real_improve / model_improve

      if model_objective > 0
        @warn "found a positive model objective increase. Debug if you see this."
        ρ = real_improve / -model_improve
      end

      if !(ρ >= solver.settings.η_2)
        tr_size = tr_size * solver.settings.t_1
      elseif ρ > solver.settings.η_3 && step_type == :boundary
        tr_size = tr_size * solver.settings.t_2
      end

      model_res = g + Jd
      model_res_norm = norm(model_res)
      real_res_norm = norm(gy)

      will_accept = 
        (ρ >= solver.settings.η_1) || 
        (ρ >= 0.0 && real_res_norm <= g_norm)

      logger(
        solver, 
        real_objective, model_objective,
        real_res_norm, model_res_norm,
        cumulative_cg_iters, tr_size, will_accept
      )
      # @assert false

      if will_accept
        # @assert false "will accept edge case"
        Uu = y
        # update_fields!(U, domain, Uu)
        g = gy
        o = energy(domain, Uu)
        g_norm = real_res_norm
        # TODO try new preconditioner upate
        happy = true
      else
        step_type = :boundary
        n_cg_iters = 0
      end

      # TODO not sure what to do here
      if n_cg_iters >= solver.settings.max_cg_iters ||
        cumulative_cg_iters >= solver.settings.max_cumulative_cg_iters
        attempt = 1 # TODO make this a global
        
        while attempt < 10
          @info "Updating preconditioner, attempt = $attempt"
          try
            P = cholesky(K; shift=10.0^(-5 + attempt))
          catch e
            @info "Failed to factor preconditioner. Attempting again"
            attempt += 1
          else
            break
          end
        end
      end
    end
  end

  # @assert false "Reached maximum number of trust region iterations"
  @warn "Reached maximum number of trust region iterations. Accepting but be careful!"
end
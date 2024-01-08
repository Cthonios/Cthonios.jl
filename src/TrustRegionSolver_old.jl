

function is_on_boundary(step_type)
  return step_type == :boundary || step_type == :negative_curvature
end

function project_to_boundary(z, d, tr_size, zz)
  dd = dot(d, d)
  zd = dot(z, d)
  tau = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
  return z + tau * d
end

function preconditioned_project_to_boundary(z, d, tr_size, zz, mult_by_approx_hessian)
  # Pd = precond * d # TODO this could be an error
  # Pd = precond \ d
  Pd = mult_by_approx_hessian(d)
  dd = dot(d, Pd)
  zd = dot(z, Pd)
  tau = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
  return z + tau * d
end

function project_to_boundary_with_coeffs(z, d, tr_size, zz, zd, dd)
  tau = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
  return z + tau * d
end

function update_step_length_squared(alpha, zz, zd, dd)
  return zz + 2 * alpha * zd + alpha * alpha * dd
end

function cg_inner_products_preconditioned(alpha, beta, zd, dd, rPr, z, d)
  zd = beta * (zd + alpha * dd)
  dd = rPr * beta * beta * dd
  return zd, dd
end

function cg_inner_products_unpreconditioned(alpha, beta, zd, dd, rPr, z, d)
  zd = dot(z, d)
  dd = dot(d, d)
  return zd, dd
end

function solve_trust_region_minimization(x, r, K, hess_vec_func, precond, tr_size, settings)
  # minimize r * z + 0.5 * z * J * z

  z  = 0.0 * x
  zz = 0.0

  cg_inexact_rel_tol = settings.cg_inexact_solve_ratio
  cg_tol_squared     = maximum((settings.cg_tol^2, cg_inexact_rel_tol^2 * dot(r, r)))
  if dot(r, r) < cg_tol_squared
    return z, z, :interior, 0
  end

  # Pr           = precond * r
  # Pr           = precond \ r
  # Pr           = K * r
  Pr           = K \ r
  d            = -Pr
  cauchy_point = d
  rPr          = dot(r, Pr)

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

    zNp1  = z + alpha * d
    zzNp1 = update_step_length_squared(alpha, zz, zd, dd)
    
    if curvature <= 0.0
      @warn "negative curvature"
      z_out = project_to_boundary_with_coeffs(z, d, tr_size, zz, zd, dd)
      return z_out, cauchy_point, :negative_curvature, i + 1
    end

    # @info "\n\n\n $zzNp1"
    if zzNp1 > tr_size^2
      @warn "projecting to boundary"
      z_out = project_to_boundary_with_coeffs(z, d, tr_size, zz, zd, dd)
      return z_out, cauchy_point, :boundary, i + 1
    end

    z = zNp1 # z + alpha * d

    r  += alpha * hess_vec_func(d)
    # Pr = precond * r
    # Pr = precond \ r
    # Pr = K * r
    Pr = K \ r
    
    rPrNp1 = dot(r, Pr)

    if dot(r, r) < cg_tol_squared
      return z, cauchy_point, :interior, i + 1
    end

    beta = rPrNp1 / rPr
    rPr = rPrNp1
    d = -Pr + beta * d

    zz = zzNp1

    zz, dd = cg_inner_products(alpha, beta, zd, dd, rPr, z, d)
  end

  return z, cauchy_point, :interior, settings.max_cg_iters
end

function dogleg_step(cp, newton_p, tr_size, mat_mul)
  cc = dot(cp, mat_mul(cp))
  nn = dot(newton_p, mat_mul(newton_p))
  tt = tr_size * tr_size

  # @show cc nn tt
  if cc >= tt # return cauchy point if it extend outside the tr
    @warn "outside trust region - returning cauchy point"
    return cp * sqrt(tt / cc)
  end

  if cc > nn # return cauchy point? seems the preconditioner was not accurate?
    @info "cp outside newton preconditioner likely inaccurate"
    return cp
  end 

  if nn > tt # on the dogleg (we have nn >= cc and tt >= cc)
    @warn "on the dogleg"
    return preconditioned_project_to_boundary(cp, newton_p - cp, tr_size, cc, mat_mul)
  end

  return newton_p
end

Base.@kwdef mutable struct TrustRegionSolverSettings <: SolverSettings
  t_1::Float64                                  = 0.25
  t_2::Float64                                  = 1.75
  η_1::Float64                                  = 1e-10
  η_2::Float64                                  = 0.1
  η_3::Float64                                  = 0.5
  max_trust_iters::Int64                        = 100
  tol::Float64                                  = 1e-8
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

function TrustRegionSolverSettings(tol::Float64)
  settings = TrustRegionSolverSettings()
  settings.tol    = tol
  settings.cg_tol = 0.2 * tol
  return settings
end

struct TrustRegionSolver{V <: AbstractVector{<:Number}}
  settings::TrustRegionSolverSettings
  ΔUu::V
end

function TrustRegionSolver(input_settings::D, domain::Domain) where D <: Dict
  settings = TrustRegionSolverSettings()
  ΔUu = zeros(Float64, length(domain.dof.unknown_indices))
  return TrustRegionSolver(settings, ΔUu)
end


function banner(
  real_obj, model_obj,
  real_res, model_res,
  cg_iters, tr_size,
  on_boundary, will_accept,
  settings
)

  @info @sprintf "O = %1.6e  O_M = %1.6e  ||R|| = %1.6e  ||R_R|| = %1.6e  CG = %3i  TR = %1.6e" real_obj model_obj real_res model_res cg_iters tr_size
end

function is_converged(
  ::TrustRegionSolver,
  x, 
  real_obj, model_obj,
  real_res, model_res,
  cg_iters, tr_size,
  settings
)

  gg = dot(real_res, real_res)
  # @show gg
  # @show settings.tol^2
  if gg < settings.tol^2
    model_res_norm = norm(model_res)
    real_res_norm  = sqrt(gg)
    banner(real_obj, model_obj, real_res_norm, model_res_norm, cg_iters, tr_size, :interior, true, settings)

    if settings.check_stability
      @assert false "not supported yet"
    end

    return true
  end

  return false
end

function solve!(
  domain::Domain, 
  solver::TrustRegionSolver,
  x::V
) where V <: AbstractVector{<:Number}

  # unpack settings
  settings          = solver.settings
  tr_size           = settings.tr_size
  tried_new_precond = false

  # initial objective and residual
  o_init = energy(domain, x)
  o = o_init
  g = residual(domain, x)
  g_norm = norm(g)
  @info "Initial objective = $o"
  @info "Initial residual  = $g_norm"

  if is_converged(solver, x, 0.0, 0.0, g, g, 0, tr_size, settings)
    return x
  end

  cumulative_cg_iters = 0

  for i in 1:settings.max_trust_iters
    if settings.use_incremental_objective
      # incremental_objective = d -> 0.5 * ((g + ))
      # @assert false "unsupported right now"
      incremental_objective = d -> 0.5 * dot(g + residual(domain, x + d), d)
    else
      incremental_objective = d -> energy(domain, x + d) - o
    end

    K = stiffness(domain, x)
    precond = cholesky(K)

    hess_vec_fun = v -> K * v

    if settings.use_preconditioned_inner_product_for_cg
      # @assert false "unsupported right now"
      @show "here in this if"
      mult_by_approx_hessian = x -> K * x # TODO maybe not right
      # mult_by_approx_hessian = x -> precond \ x
    else
      mult_by_approx_hessian = x -> x
    end

    # now the meat

    gKg = dot(g, hess_vec_fun(g))

    if gKg > 0
      alpha = -dot(g, g) / gKg
      cauchy_point = alpha * g
      cauchy_point_norm_squared = dot(cauchy_point, mult_by_approx_hessian(cauchy_point))
    else
      cauchy_point = -g * (tr_size / sqrt(dot(g, mult_by_approx_hessian(g))))
      cauchy_point_norm_squared = tr_size * tr_size
      @info "negative curavture unpreconditioned cauchy point direction found."
    end

    if cauchy_point_norm_squared >= tr_size * tr_size
      @info "unpreconditioned gradient cauchy point outside trust region at dist = $(sqrt(cauchy_point_norm_squared))."
      cauchy_point = (tr_size / sqrt(cauchy_point_norm_squared)) * cauchy_point
      cauchy_point_norm_squared = tr_size * tr_size
      q_newton_point = cauchy_point
      step_type = :boundary
      cg_iters = 1
    else
      q_newton_point, _, step_type, cg_iters = solve_trust_region_minimization(x, g, K, hess_vec_fun, precond, tr_size, settings)
    end

    cumulative_cg_iters += cg_iters

    tr_size_used = tr_size

    happy_about_tr_size = false

    while !happy_about_tr_size
      # take dogleg step
      @info "taking dogleg step"
      d = dogleg_step(cauchy_point, q_newton_point, tr_size, mult_by_approx_hessian)

      Jd = hess_vec_fun(d)
      dJd = dot(d, Jd)
      # model_objective = dot(g, d .+ 0.5 * dJd)
      model_objective = dot(g, d) + 0.5 * dJd

      y = x + d
      real_objective = incremental_objective(d)
      gy = residual(domain, y)

      if is_converged(solver, y, real_objective, model_objective, gy, g + Jd, cg_iters, tr_size_used, settings)
        @info "converged here"
        return y
      end

      model_improve = -model_objective
      real_improve  = -real_objective

      rho = real_improve / model_improve

      if model_objective > 0.0
        @info "Found a positive model objective increase. Debug if you see this"
        rho = real_improve / -model_improve
      end

      @info "rho = $rho"
      @info "$(rho >= settings.η_2)"
      @info "step type = $step_type"
      @info "$(is_on_boundary(step_type))"

      if !(rho >= settings.η_2) # write it this way to handle NaNs
        tr_size *= settings.t_1
      elseif rho > settings.η_3 && is_on_boundary(step_type)
        tr_size *= settings.t_2
      end      

      model_res = g + Jd
      model_res_norm = norm(model_res)
      real_res_norm  = norm(gy)

      @show real_res_norm
      @show g_norm

      will_accept = rho >= settings.η_1 || (rho >= -0.0 && real_res_norm <= g_norm)

      @show will_accept
      banner(real_objective, model_objective, real_res_norm, model_res_norm, cg_iters, tr_size_used, step_type, will_accept, settings)
      
      if will_accept
        @info "accepting"
        x = y
        g = gy
        o = energy(domain, x)
        incremental_objective = d -> energy(domain, x + d) - o
        g_norm = real_res_norm
        tried_new_precond = false
        @show happy_about_tr_size = true
      else
        step_type = :boundary
        cg_iters = 0
      end

      if cg_iters >= settings.max_cg_iters || cumulative_cg_iters >= settings.max_cumulative_cg_iters
        # update preconditioner here TODO
        # @assert false
        @info "Updating preconditioner"
        attempt = 1
        # K = domain.assembler.K[domain.dof.is_unknown, domain.dof.is_unknown]
        dAbs = abs.(diag(K))
        shift = 10^(-5. + attempt)
        # precond = LinearMap(K .+ shift* dAbs * I)
        # precond = cholesky(K; shift=shift * dAbs)
        # temp = K + shift * dAbs .* I
        temp = K
        for n in size(K, 1)
          temp[n, n] += shift * dAbs[n]
        end
        precond = cholesky(temp)
        cumulative_cg_iters = 0
      end

      tr_size_used = tr_size

      if tr_size < settings.min_tr_size
        if !tried_new_precond
          @info "The trust region is too small, updating precond and trying again."
          dAbs = abs.(diag(K))
          shift = 10^(-5. + attempt)
          # precond = LinearMap(K .+ shift* dAbs * I)
          precond = cholesky(K; shift=shift * dAbs)
          cumulative_cg_iters = 0
          tried_new_precond = true
          happy_about_tr_size = true
          tr_size = settings.tr_size
        else
          @info "The trust region is still too small. Accepting, but be careful."
          return x
        end
      end
      # @assert false "made it here"
    end
  end

  @info "Reached the maximum number of trust region iterations."
  if settings.check_stability
    # add stability check here
  end

  return x
end

function solve!(
  domain::Domain, 
  solver::TrustRegionSolver;
  use_warm_start::Bool=false,
  update_precond=true
)
  scaling     = 1.0 # TODO maybe change this?
  inv_scaling = 1.0

  x0 = solver.ΔUu

  FiniteElementContainers.update_fields!(domain.U, domain.dof, x0)

  x_bar0 = scaling * x0

  if use_warm_start
    @assert false "warm start not implemented yet"
  else
    # do nothing here
  end

  if update_precond
    @info "Updating preconditioner"
    K = stiffness(domain, x_bar0)
    precond = cholesky(K)
  end

  x_bar = solve!(domain, solver, x_bar0)
  solver.ΔUu .= x_bar
end
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
  F <: NodalField,
  S <: TrustRegionSolverSettings
} <: AbstractNonlinearSolver{L, O, U}
  preconditioner::L
  objective::O
  ΔUu::U
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
  settings::S
end

"""
$(TYPEDSIGNATURES)
TODO figure out which scratch arrays can be nixed
"""
function TrustRegionSolver(objective, p; use_warm_start=true)
  domain = objective.domain
  settings       = TrustRegionSolverSettings() # TODO add non-defaults
  # TODO eventually write a custom linear solver for this one
  preconditioner = CholeskyPreconditioner(objective, p)
  ΔUu            = create_unknowns(domain)
  o              = zeros(1)
  g              = create_fields(domain) # gradient
  Hv             = create_fields(domain) # hessian-vector product
  cauchy_point   = create_unknowns(domain)
  q_newton_point = create_unknowns(domain)
  d              = create_unknowns(domain)
  y_scratch_1    = create_unknowns(domain)
  y_scratch_2    = create_unknowns(domain)
  y_scratch_3    = create_unknowns(domain)
  y_scratch_4    = create_unknowns(domain)
  # TODO remove this hardcoded thing
  use_warm_start = false
  return TrustRegionSolver(
    preconditioner, objective,
    ΔUu,
    o, g, Hv,
    cauchy_point, q_newton_point, d,
    y_scratch_1, y_scratch_2, y_scratch_3, y_scratch_4,
    use_warm_start,
    settings
  )
end

function Base.show(io::IO, solver::TrustRegionSolver)
  print(io, 
    "    TrustRegionSolver\n", 
    "    $(solver.settings)\n",
    "    $(solver.linear_solver)\n"
  )
end

function logger(
  ::TrustRegionSolver,
  real_obj, model_obj,
  real_res, model_res,
  cg_iters, tr_size, accept, step_type
)

  str = @sprintf "O = %1.4e  O_M = %1.4e  ||R|| = %1.4e  ||R_R|| = %1.4e  CG = %3i  TR = %1.4e accept = %s %s" real_obj model_obj real_res model_res cg_iters tr_size accept step_type
  @info str
end

function is_converged(
  solver::TrustRegionSolver,
  real_obj, model_obj,
  real_res, model_res,
  cg_iters, tr_size, step_type
)

  gg = dot(real_res, real_res)

  if gg < solver.settings.tol^2
    model_res_norm = norm(model_res)
    real_res_norm  = sqrt(gg)
    logger(solver, real_obj, model_obj, real_res_norm, model_res_norm, cg_iters, tr_size, true, step_type)

    if solver.settings.check_stability
      @assert false "not supported yet"
    end

    return true
  end

  return false
end

"""
$(TYPEDSIGNATURES)
"""
function calculate_cauchy_point!(cauchy_point, solver, P, g, Hg, tr_size)
  y_scratch_1 = solver.y_scratch_1
  gHg = dot(g, Hg)
  if gHg > 0
    α = -dot(g, g) / gHg
    @. cauchy_point = α * g
    ldiv!(y_scratch_1, P, cauchy_point)
    cauchy_point_norm_squared = dot(cauchy_point, y_scratch_1)
  else
    ldiv!(y_scratch_1, P, g)
    @. cauchy_point = -g * (tr_size / sqrt(dot(g, y_scratch_1)))
    @info "negative curavture unpreconditioned cauchy point direction found."
  end
  return cauchy_point_norm_squared
end

"""
$(TYPEDSIGNATURES)
"""
function update_tr_size(solver, model_objective, real_objective, step_type, tr_size, real_res_norm, g_norm)
  ρ = -model_objective / -real_objective

  if model_objective > 0
    @warn "found a positive model objective increase. Debug if you see this."
    ρ = -real_objective / model_objective
  end

  if !(ρ >= solver.settings.η_2)
    tr_size = tr_size * solver.settings.t_1
  elseif ρ > solver.settings.η_3 && step_type == :boundary
    tr_size = tr_size * solver.settings.t_2
  end

  will_accept = 
    (ρ >= solver.settings.η_1) || 
    (ρ >= 0.0 && real_res_norm <= g_norm)

  return tr_size, will_accept
end

"""
$(TYPEDSIGNATURES)
minimize r * z + 0.5 * z * J * z
"""
function minimize_trust_region_sub_problem!(
  z,
  solver::TrustRegionSolver, P,#domain, #common::CthoniosCommon,
  x::V, p, r_in, #cauchy_point, # cauhcy point is really a scratch array to reduce an allocation
  tr_size::Float64
) where V <: AbstractVector

  Pr   = solver.y_scratch_1
  d    = solver.y_scratch_2
  r    = solver.y_scratch_3
  zNp1 = solver.y_scratch_4

  # zero some stuff out and make other scratch vars
  z .= zero(eltype(z))
  zz = zero(eltype(x))
  zd = zero(eltype(x))
  r .= r_in

  cg_tol_squared = max(
    solver.settings.cg_inexact_solve_ratio^2 * dot(r, r),
    solver.settings.cg_tol^2
  )

  ldiv!(Pr, P, r)

  @. d = -Pr

  rPr = dot(r, Pr)

  dd = rPr

  for i in 1:solver.settings.max_cg_iters
    Hd = hvp(solver, x, p, d)

    curvature = dot(d, Hd)
    α = rPr / curvature

    @. zNp1 = z + α * d
    zzNp1 = zz + 2.0 * α * zd + α * α * dd

    # TODO
    if curvature <= zero(eltype(x))
      τ = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
      @. z = z + τ * d
      return :negative_curvature, i
    end

    # TODO
    if zzNp1 > tr_size^2
      # @assert false "other edge case not handled"
      # TODO projec to boundary
      τ = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
      @. z = z + τ * d
      return :boundary, i
    end

    @. z = zNp1

    @. r = r + α * Hd

    ldiv!(Pr, P, r)

    rPrNp1 = dot(r, Pr)
    
    if dot(r, r) < cg_tol_squared
      return :interior, i
    end

    β = rPrNp1 / rPr
    rPr = rPrNp1
    @. d = -Pr + β * d

    zz = zzNp1
    zd = β * (zd + α * dd)
    dd = rPr * β * β * dd
  end
  # end
  return :interior, solver.settings.max_cg_iters
end

"""
$(TYPEDSIGNATURES)
"""
function dog_leg_step!(d, solver, P, cauchy_point, q_newton_point, tr_size)
  y_scratch_1 = solver.y_scratch_1
  y_scratch_2 = solver.y_scratch_2

  ldiv!(y_scratch_1, P, cauchy_point)
  cc = dot(cauchy_point, y_scratch_1)
  ldiv!(y_scratch_1, P, q_newton_point)
  nn = dot(q_newton_point, y_scratch_1)

  tt = tr_size * tr_size

  # return cauchy point if it extends outside the tr
  if cc >= tt
    d .= sqrt(tt / cc) * cauchy_point
    return nothing
  end

  # return cauchy point? seems the preconditioner was not accurate
  if cc > nn
    @warn "cp outside newton, preconditioner likely inaccurate"
    d .= cauchy_point
    return nothing
  end

  # on the dogleg
  if nn > tt
    @. y_scratch_2 = q_newton_point - cauchy_point
    ldiv!(y_scratch_1, P, y_scratch_2)
    dd = dot(y_scratch_2, y_scratch_1)
    zd = dot(cauchy_point, y_scratch_1)

    τ  = (sqrt((tr_size^2 - cc) * dd + zd^2) - zd) / dd
    @. d = cauchy_point + τ * (y_scratch_2)
    return nothing
  end

  d .= q_newton_point
  return nothing
end

"""
$(TYPEDSIGNATURES)
TODO
Eventually map this to the interface
"""
function solve!(solver::TrustRegionSolver, Uu, p)
  @info "TrustRegionSolver"
  # unpack some stuff
  ΔUu            = solver.ΔUu
  cauchy_point   = solver.cauchy_point
  q_newton_point = solver.q_newton_point
  d              = solver.d

  # unpack some solver settings
  tr_size = solver.settings.tr_size

  # saving initial objective and gradient
  o = objective(solver, Uu, p)
  g = gradient(solver, Uu, p)
  g_norm = norm(g)

  @info @sprintf "Initial objective = %1.6e" o[1]
  @info @sprintf "Initial residual  = %1.6e" g_norm

  # preconditioner
  P = solver.preconditioner
  update_preconditioner!(P, solver.objective, Uu, p)

  # begin loop over tr iters
  cumulative_cg_iters = 0
  tried_new_precond = false

  for n in 1:solver.settings.max_trust_iters
    @info "Trust region iter = $n"
    # need hvp
    Hg = hvp(solver, Uu, p, g)

    # check for negative curvature
    cauchy_point_norm_squared = calculate_cauchy_point!(cauchy_point, solver, P, g, Hg, tr_size)

    # check if outside trust region
    if cauchy_point_norm_squared >= tr_size * tr_size
      @info "unpreconditioned gradient cauchy point outside trust region at dist = $(sqrt(cauchy_point_norm_squared))."
      @. cauchy_point = (tr_size / sqrt(cauchy_point_norm_squared)) * cauchy_point
      cauchy_point_norm_squared = tr_size * tr_size
      q_newton_point .= cauchy_point
      step_type = :boundary
      n_cg_iters = 1
    else
      step_type, n_cg_iters = 
        minimize_trust_region_sub_problem!(q_newton_point, solver, P, Uu, p, g, tr_size)
    end

    cumulative_cg_iters += n_cg_iters

    happy = false
    while !happy
      dog_leg_step!(d, solver, P, cauchy_point, q_newton_point, tr_size)

      Jd = hvp(solver, Uu, p, d)
      dJd = dot(d, Jd)

      model_objective = dot(g, d) + 0.5 * dJd

      y = solver.y_scratch_1
      @. y = Uu + d

      real_objective = Cthonios.objective(solver, y, p)
      gy = Cthonios.gradient(solver, y, p)
      real_objective = real_objective - o

      if is_converged(solver, real_objective, model_objective, gy, g + Jd, n_cg_iters, tr_size, step_type)
        @info "Converged"
        @. ΔUu = Uu - solver.y_scratch_1
        Uu .= solver.y_scratch_1
        return nothing
      end

      model_res = solver.y_scratch_2
      @. model_res = g + Jd
      model_res_norm = norm(model_res)
      real_res_norm = norm(gy)

      tr_size, will_accept = update_tr_size(
        solver, 
        model_objective, real_objective, 
        step_type, tr_size, 
        real_res_norm, g_norm
      )

      logger(
        solver, 
        real_objective, model_objective,
        real_res_norm, model_res_norm,
        cumulative_cg_iters, tr_size, will_accept, step_type
      )

      if will_accept
        Uu .= y
        g .= gy
        o = objective(solver, Uu, p)
        g_norm = norm(g)
        happy = true
      else
        step_type = :boundary
        n_cg_iters = 0
      end

      # TODO not sure what to do here
      if n_cg_iters >= solver.settings.max_cg_iters ||
         cumulative_cg_iters >= solver.settings.max_cumulative_cg_iters
        update_preconditioner!(P, solver.objective, Uu, p)
      end

      # TODO still need to do some other stuff with updating preconditioner
      if tr_size < solver.settings.min_tr_size
        if !tried_new_precond
          @info "Updating preconditioner due to small trust region size"
          update_preconditioner!(P, solver.objective, Uu, p)
          cumulative_cg_iters = 0
          tried_new_precond = true
          happy = true
          tr_size = solver.settings.tr_size
        else
          @info "The trust region is too small. Accepting but be careful."
          return nothing
        end
      end
    end
  end
  @assert false "Reached maximum number of trust region iterations."

end

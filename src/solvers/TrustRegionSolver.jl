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
    print(io, "        ", rpad(name, 39), " = ", getfield(settings, name), "\n")
  end
  print(io, "\n")
end

struct TrustRegionSolver{
  S <: TrustRegionSolverSettings,
  L <: AbstractLinearSolver,
  V <: AbstractVector,
  F <: NodalField,
  C1, C2, C3, C4, C5
} <: NonlinearSolver{S, L, V}
  settings::S
  linear_solver::L
  Uu::V
  ΔUu::V
  # o::V
  # g::V
  Hv::V
  Hv_field::F
  #
  cauchy_point::C1
  y::C2
  d::C3
  z::C4
  r::C5
end

function TrustRegionSolver(input_settings::D, domain::QuasiStaticDomain) where D <: Dict{Symbol, Any}
  settings      = TrustRegionSolverSettings() # TODO add non-defaults
  linear_solver = setup_linear_solver(input_settings[Symbol("linear solver")], domain)
  Uu            = create_unknowns(domain)
  ΔUu           = create_unknowns(domain)
  # o             = zeros(Float64, 1)
  # g             = create_unknowns(domain)
  Hv            = create_unknowns(domain)
  Hv_field      = create_fields(domain)
  # 
  # cache arrays for sub problem
  #
  cauchy_point  = create_unknowns(domain)
  y             = create_unknowns(domain)
  d             = create_unknowns(domain)
  z             = create_unknowns(domain)
  r             = create_unknowns(domain)
  return TrustRegionSolver(settings, linear_solver, Uu, ΔUu, Hv, Hv_field, cauchy_point, y, d, z, r)
end

function Base.show(io::IO, solver::TrustRegionSolver)
  print(io, "    TrustRegionSolver\n", 
        "    $(solver.settings)\n",
        "    $(solver.linear_solver)\n")
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

function objective(solver, domain, common, u)
  @timeit timer(common) "Energy" begin
    energy!(solver.linear_solver, domain, u)
    o = domain.domain_cache.Π[1]
  end
  return o
end

function objective_and_grad(solver, domain, common, u)
  @timeit timer(common) "Energy and internal force" begin
    energy_and_internal_force!(solver.linear_solver, domain, u)
    o = domain.domain_cache.Π[1]
    # g = @views solver.linear_solver.assembler.residuals[domain.dof.unknown_dofs]
    g = @views domain.domain_cache.f[domain.dof.unknown_dofs]
  end
  return o, g
end

function hvp(solver, domain, common, u, v)
  @timeit timer(common) "Stiffness action" begin
    stiffness_action!(solver.Hv_field, domain, u, v)
    Hv = @views solver.Hv_field[domain.dof.unknown_dofs]
  end
  return Hv
end

function hvp!(Hv, solver, domain, common, u, v)
  @timeit timer(common) "Stiffness action" begin
    stiffness_action!(solver.Hv_field, domain, u, v)
    Hv .= @views solver.Hv_field[domain.dof.unknown_dofs]
  end
  return nothing
end

function hessian(solver, domain, common, u)
  @timeit timer(common) "Hessian" begin
    stiffness!(solver, domain, u)
    K = SparseArrays.sparse!(solver.linear_solver.assembler) #|> Symmetric
  end
  return K
end

function preconditioner(solver, domain, common, Uu)
  @timeit timer(common) "Preconditioner" begin
    H = hessian(solver, domain, common, Uu)
    # P = ldl(H)
    ldl_factorize!(H, solver.linear_solver.factor)
    return P
    # attempt = 1 # TODO make this a global
    # while attempt < 10
    #   @info "Updating preconditioner, attempt = $attempt"
    #   try
    #     # P = cholesky(H; shift=10.0^(-5 + attempt))
    #     P = ldl(H)
    #     # P = opCholesky(H)
    #     # P = CholeskyPreconditioner(H)
    #     # P = DiagonalPreconditioner(H)
    #     # NU = length(domain.dof.unknown_dofs)
    #     # op = LinearOperator(
    #     #   Float64, NU, NU, true, true,
    #     #   (Hv, v) -> hvp!(Hv, solver, domain, common, Uu, v)
    #     # )
    #     # P = InvPreconditioner(op)
    #     return P
    #   catch e
    #     @info e
    #     @info "Failed to factor preconditioner. Attempting again"
    #     attempt += 1
    #   else
    #     break
    #   end
    # end
    # @assert false "Failed to factorize hessian after 10 attempts"
  end
end

function preconditioner!(solver, domain, common, Uu)
  @timeit timer(common) "Preconditioner" begin
    H = hessian(solver, domain, common, Uu)
    ldl_factorize!(H, solver.linear_solver.factorization)
  end
  return nothing
end

function calculate_cauchy_point!(solver, common, g, Hg, tr_size)
  @timeit timer(common) "Cauchy point" begin
    P = solver.linear_solver.factorization
    gHg = dot(g, Hg)
    if gHg > 0
      α = -dot(g, g) / gHg
      mul!(solver.cauchy_point, α, g)
      ldiv!(solver.y, P, solver.cauchy_point)
      cauchy_point_norm_squared = dot(solver.cauchy_point, solver.y)
    else
      ldiv!(solver.y, P, g)
      mul!(solver.cauchy_point, -(tr_size / sqrt(dot(g, solver.y))), g)
      cauchy_point_norm_squared = tr_size * tr_size
      @info "negative curavture unpreconditioned cauchy point direction found."
    end
  end
  return cauchy_point_norm_squared
end

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
minimize r * z + 0.5 * z * J * z
"""
function minimize_trust_region_sub_problem(
  solver::TrustRegionSolver, domain, common::CthoniosCommon,
  x::V, r_in::V, P, 
  # scratch arrays
  # TODO wrap these in a cache struct
  # y_scratch, d, z, r, # cauhcy point is really a scratch array to reduce an allocation
  tr_size::Float64
) where V <: AbstractVector

  @timeit timer(common) "Subproblem" begin

    solver.z .= zero(eltype(solver.z))
    solver.r .= r_in

    cg_tol_squared = max(
      solver.settings.cg_inexact_solve_ratio^2 * dot(solver.r, solver.r),
      solver.settings.cg_tol^2
    )

    @timeit timer(common) "Multiply by approximate hessian" begin
      ldiv!(solver.y, P, solver.r)
    end

    # d = -Pr
    @. solver.d = -solver.y
    # cauchy_point = zeros(eltype(d), length(d))
    # cauchy_point = copy(d)
    cauchy_point = solver.d

    rPr = dot(solver.r, solver.y)

    zz = zero(eltype(x))
    zd = zero(eltype(x))

    dd = rPr

    for i in 1:solver.settings.max_cg_iters
      Hd = hvp(solver, domain, common, x, solver.d)

      # curvature = dot(d, K * d)
      curvature = dot(solver.d, Hd)
      α = rPr / curvature

      # @. zNp1 = z + α * d
      @. solver.y = solver.z + α * solver.d
      zzNp1 = zz + 2.0 * α * zd + α * α * dd

      # TODO
      if curvature <= zero(eltype(x))
        τ = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
        return solver.z + τ * solver.d, cauchy_point, :negative_curvature, i
      end

      # TODO
      if zzNp1 > tr_size^2
        # @assert false "other edge case not handled"
        # TODO projec to boundary
        τ = (sqrt((tr_size^2 - zz) * dd + zd^2) - zd) / dd
        return solver.z + τ * solver.d, cauchy_point, :boundary, i
      end

      # z = zNp1
      @. solver.z = solver.y

      @. solver.r = solver.r + α * Hd

      @timeit timer(common) "Multiply by approximate hessian" begin
        ldiv!(solver.y, P, solver.r)
      end

      rPrNp1 = dot(solver.r, solver.y)
      
      if dot(solver.r, solver.r) < cg_tol_squared
        return solver.z, cauchy_point, :interior, i
      end

      β = rPrNp1 / rPr
      rPr = rPrNp1
      @. solver.d = -solver.y + β * solver.d

      zz = zzNp1
      zd = β * (zd + α * dd)
      dd = rPr * β * β * dd
    end
  end
  return solver.z, cauchy_point, :interior, solver.settings.max_cg_iters
end

function dog_leg_step(solver, common, q_newton_point, tr_size, P)
  @timeit timer(common) "Doglog step" begin
    # old
    # cc = dot(cauchy_point, P \ cauchy_point)
    # nn = dot(q_newton_point, P \ q_newton_point)
    ldiv!(solver.y, P, solver.cauchy_point)
    cc = dot(solver.cauchy_point, solver.y)
    ldiv!(solver.y, P, q_newton_point)
    nn = dot(q_newton_point, solver.y)

    tt = tr_size * tr_size

    # return cauchy point if it extends outside the tr
    if cc >= tt
      return solver.cauchy_point * sqrt(tt / cc)
    end

    # return cauchy point? seems the preconditioner was not accurate
    if cc > nn
      @warn "cp outside newton, preconditioner likely inaccurate"
      return solver.cauchy_point
    end

    # on the dogleg
    if nn > tt
      mul!(solver.y, P, q_newton_point - solver.cauchy_point)      
      dd = dot(q_newton_point - solver.cauchy_point, solver.y)
      zd = dot(solver.cauchy_point, solver.y)
      τ  = (sqrt((tr_size^2 - cc) * dd + zd^2) - zd) / dd
      return solver.cauchy_point + τ * (q_newton_point - solver.cauchy_point)
    end
  end
  return q_newton_point
end

function solve!(
  solver::TrustRegionSolver, domain::QuasiStaticDomain,
  common::CthoniosCommon
)

  # unpack cached arrays from solver and domain
  Uu, ΔUu = solver.Uu, solver.ΔUu

  # unpack some solver settings
  tr_size = solver.settings.tr_size

  # calculate initial residual and objective
  @timeit timer(common) "Energy, internal force, and stiffness" begin
    energy_internal_force_and_stiffness!(solver.linear_solver, domain, Uu)
  end

  # saving initial objective and gradient
  o = domain.domain_cache.Π[1]
  g = domain.domain_cache.f[domain.dof.unknown_dofs]
  g_norm = norm(g)
  y = create_unknowns(domain)
  model_res = similar(g)
  model_res .= zero(eltype(g))

  # preconditioner
  preconditioner!(solver, domain, common, Uu)
  P = solver.linear_solver.factorization

  @info @sprintf "Initial objective = %1.6e" o
  @info @sprintf "Initial residual  = %1.6e" g_norm

  # begin loop over tr iters
  cumulative_cg_iters = 0
  for n in 1:solver.settings.max_trust_iters

    # check for convergence at beginning of iteration
    if is_converged(solver, Uu, 0.0, 0.0, g, g, 0, tr_size)
      @info "Converged - figure out what to do here"
    end

    # need hvp
    Hg = hvp(solver, domain, common, Uu, g)

    # check for negative curvature
    cauchy_point_norm_squared = calculate_cauchy_point!(solver, common, g, Hg, tr_size)

    # check if outside trust region
    if cauchy_point_norm_squared >= tr_size * tr_size
      @info "unpreconditioned gradient cauchy point outside trust region at dist = $(sqrt(cauchy_point_norm_squared))."
      cauchy_point .= (tr_size / sqrt(cauchy_point_norm_squared)) * solver.cauchy_point
      cauchy_point_norm_squared = tr_size * tr_size
      q_newton_point = solver.cauchy_point
      step_type = :boundary
      n_cg_iters = 1
    else
      q_newton_point, _, step_type, n_cg_iters = 
        minimize_trust_region_sub_problem(
          solver, domain, common, Uu, g, P, 
          tr_size
        )

      step_type = :cg
    end

    cumulative_cg_iters += n_cg_iters

    happy = false

    while !happy
      d = dog_leg_step(solver, common, q_newton_point, tr_size, P)

      Jd = hvp(solver, domain, common, Uu, d)
      dJd = dot(d, Jd)

      model_objective = dot(g, d) + 0.5 * dJd

      # y = Uu + d
      @. y = Uu + d

      real_objective, gy = objective_and_grad(solver, domain, common, y)
      real_objective = real_objective - o

      if is_converged(solver, y, real_objective, model_objective, gy, g + Jd, n_cg_iters, tr_size)
        @info "Converged"
        ΔUu .= Uu - y
        Uu .= y
        return nothing
      end

      # model_res = g + Jd
      @. model_res = g + Jd
      model_res_norm = norm(model_res)
      real_res_norm = norm(gy)

      tr_size, will_accept = update_tr_size(solver, model_objective, real_objective, step_type, tr_size, real_res_norm, g_norm)

      logger(
        solver, 
        real_objective, model_objective,
        real_res_norm, model_res_norm,
        cumulative_cg_iters, tr_size, will_accept
      )

      if will_accept
        Uu .= y
        g .= gy
        o = objective(solver, domain, common, Uu)
        g_norm = norm(g)
        happy = true
      else
        step_type = :boundary
        n_cg_iters = 0
      end

      # TODO not sure what to do here
      if n_cg_iters >= solver.settings.max_cg_iters ||
        cumulative_cg_iters >= solver.settings.max_cumulative_cg_iters
        # P = preconditioner(solver, domain, common, Uu)
        preconditioner!(solver, domain, common, Uu)
      end
    end
  end

  @assert false
end
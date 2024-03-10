function seed!(::ReverseMode, cache::Cthonios.QuasiStaticDomainCache)
  cache.X .= 0.0
  cache.Uu .= 0.0
  cache.ΔUu .= 0.0
  cache.U .= 0.0
  cache.props .= 0.0
  cache.state_old .= 0.0
  cache.state_new .= 0.0
  cache.Π .= 0.0
  cache.Πs .= 0.0
  cache.f .= 0.0
  cache.V .= 0.0
  cache.Hv .= 0.0
  return nothing
end

function warm_start!(solver, domain::QuasiStaticDomain, common, Uu)
  @info "Warm start"
  # unpack stuff
  cache = domain.domain_cache
  Δt = cache.time.Δt
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, f = cache

  # update fields in a round about way
  Uu_old = copy(Uu)
  U_old  = copy(U)
  update_fields!(U_old, domain, X, Uu_old)

  # now update time
  step!(domain)
  update_fields!(U, domain, X, Uu)

  dU = similar(U)
  dU .= 0.0
  @. dU = U_old - U

  K = hessian(solver, domain, common, Uu)

  @timeit timer(common) "JacVec" begin

    Jv = JacVec(
      z -> internal_force(sections, Δt, X, z, props, state_old, backend(common))[1],
      U
    )
  end

  @timeit timer(common) "Jv*dU" mul!(cache.Hv, Jv, dU)

  @views b = cache.Hv[domain.dof.unknown_dofs]
  @timeit timer(common) "solve" temp_Uu = K \ b

  @. cache.Uu += temp_Uu
end

function warm_start_old!(solver, domain::QuasiStaticDomain, common, Uu)
  @info "Warm start"
  # unpack stuff
  cache = domain.domain_cache
  Δt = cache.time.Δt
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, f = cache

  # update fields in a round about way
  Uu_old = copy(Uu)
  U_old  = copy(U)
  update_fields!(U_old, domain, X, Uu_old)

  # now update time
  step!(domain)
  update_fields!(U, domain, X, Uu)

  dU = similar(U)
  dU .= 0.0
  @. dU = U_old - U

  # AD setup
  dcache = similar(cache)
  seed!(Reverse, dcache)
  dcache.U .= 1.0

  # using stiffness matrix for now, switch to linear operator
  stiffness!(solver, domain, cache, Uu, backend(common))
  K = SparseArrays.sparse!(solver.linear_solver.assembler)

  autodiff(
    Forward, internal_force!,
    Duplicated(cache.f, dcache.f),
    Duplicated(cache.state_new, dcache.state_new),
    Const(sections),
    Const(Δt),
    Const(cache.X),
    Duplicated(cache.U, dU),
    Const(cache.props),
    Const(cache.state_old),
    Const(backend(common))
  )

  @views b = dcache.f[domain.dof.unknown_dofs]
  temp_Uu = K \ b

  @. cache.Uu += temp_Uu

  return nothing
end

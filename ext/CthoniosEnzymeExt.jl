module CthoniosEnzymeExt

using Cthonios
using Enzyme
using EnzymeCore
using LinearAlgebra


# method for full gradients of everything
function Cthonios.gradient(domain::Cthonios.QuasiStaticDomain, backend)
  # primals
  Uu = domain.domain_cache.Uu
  U = domain.domain_cache.U
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords
  Π = zeros(eltype(U), 1)

  # gradients
  dUu = similar(Uu)
  dU = similar(U)
  dstate = similar(state)
  dprops = similar(props)
  dX = similar(X)
  dΠ = similar(Π)

  # seeding
  dUu .= zero(eltype(dUu))
  dU .= zero(eltype(dU))
  dstate .= zero(eltype(dstate))
  dprops .= zero(eltype(dprops))
  dX .= zero(eltype(dX))
  dΠ .= one(eltype(dΠ))

  # calculate gradient
  Cthonios.gradient!(
    dΠ, dU, dUu, dstate, dprops, dX,
    Π, U, domain, Uu, state, props, X,
    backend
  )
end

function Cthonios.gradient!(
  dΠ, dU, dUu, dstate, dprops, dX,
  Π, U, domain, Uu, state, props, X,
  backend
)
  autodiff_deferred(
    Reverse, Cthonios.energy!,
    Duplicated(Π, dΠ),
    Duplicated(U, dU),
    Const(domain),
    Duplicated(Uu, dUu),
    Duplicated(state, dstate),
    Duplicated(props, dprops),
    Duplicated(X, dX),
    Const(backend)
  )
  return nothing
end

function Cthonios.residual(domain::Cthonios.QuasiStaticDomain)
  # primals
  Uu = domain.domain_cache.Uu
  U = domain.domain_cache.U
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords
  Π = zeros(eltype(U), 1)

  # gradients
  dUu = similar(Uu)
  dU = similar(U)
  dΠ = similar(Π)

  # seeding
  dUu .= zero(eltype(dUu))
  dU .= zero(eltype(dU))
  dΠ .= one(eltype(dΠ))

  # use in place method
  Cthonios.residual!(dΠ, dU, dUu, Π, U, domain, Uu, state, props, X)

  return dUu
end

function Cthonios.residual!(
  dΠ, dU, dUu,
  Π, U, domain, Uu, state, props, X
)
  autodiff(
    Reverse, Cthonios.energy!,
    Duplicated(Π, dΠ),
    Duplicated(U, dU),
    Const(domain),
    Duplicated(Uu, dUu),
    Const(state),
    Const(props),
    Const(X)
  )
  return nothing
end

function Cthonios.residual_dot_v!(
  Rv, dΠ, dU, dUu,
  Π, U, domain, Uu, state, props, X, v
)
  Cthonios.residual!(dΠ, dU, dUu, Π, U, domain, Uu, state, props, X)
  Rv[1] = dot(dUu, v)
  return nothing
end

function Cthonios.hvp(domain::Cthonios.QuasiStaticDomain, v)
  # primals
  Uu = domain.domain_cache.Uu
  U = domain.domain_cache.U
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords
  Π = zeros(eltype(U), 1)
  Rv = zeros(eltype(U), 1)

  # gradients
  bUu = similar(Uu)
  bU = similar(U)
  bΠ = similar(Π)
  bRv = similar(Rv)

  dUu = similar(Uu)
  dU = similar(U)
  dΠ = similar(Π)
  dRv = similar(Rv)

  dbUu = similar(Uu)
  dbU = similar(U)
  dbΠ = similar(Π)
  dbRv = similar(Rv)

  # seeding
  bUu .= zero(eltype(bUu))
  bU .= zero(eltype(bU))
  bΠ .= one(eltype(bΠ))
  bRv .= one(eltype(bRv))

  dUu .= one(eltype(dUu))
  dU  .= zero(eltype(dU))
  dΠ .= zero(eltype(dΠ))
  dRv .= zero(eltype(dRv))

  dbUu .= zero(eltype(dbUu))
  dbU .= zero(eltype(dbU))
  dbΠ .= zero(eltype(dbΠ))
  dbRv .= zero(eltype(dbRv))

  # use in place method
  # Cthonios.residual_dot_v!(Rv, bΠ, bU, bUu, Π, U, domain, Uu, state, props, X, v)
  # Rv

  # autodiff(
  #   Forward,
  #   Cthonios.residual_dot_v!,
  #   Duplicated(Duplicated(Rv, bRv), Duplicated(dRv, dbRv)),
  #   Duplicated(Duplicated(Π, bΠ), Duplicated(dΠ, dbΠ)),
  # )

  # autodiff(
  #   Forward,
  #   (a, b, c, d, e, f, g) -> autodiff_deferred(Reverse, Cthonios.energy!, a, b, c, d, e, f, g),
  #   Duplicated(Duplicated(Π, bΠ), Duplicated(dΠ, dbΠ)),
  #   Duplicated(Duplicated(U, bU), Duplicated(dU, dbU)),
  #   Const(domain),
  #   Duplicated(Duplicated(Uu, bUu), Duplicated(dUu, dbUu)),
  #   Const(state),
  #   Const(props),
  #   Const(X)  
  # )

  autodiff(
    Forward,
    Cthonios.energy!,
    Duplicated(Π, dΠ),
    Duplicated(U, dU),
    Const(domain),
    Duplicated(Uu, dUu),
    Const(state),
    Const(props),
    Const(X)  
  )
  dU
end

end # module
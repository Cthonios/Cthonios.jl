# lazy method for developers in the REPL
function energy(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return energy(domain, Uu, state, props, X)
end    

function internal_force(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return internal_force(domain, Uu, state, props, X)
end   

function stiffness(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return stiffness(domain, Uu, state, props, X)
end

function internal_force_and_stiffness(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return internal_force_and_stiffness(domain, Uu, state, props, X)
end

# out of place methods with differentiable stuff exposed
function energy(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  Π = zeros(eltype(Uu), 1)
  energy!(Π, U, domain, Uu, state, props, X)
  return Π
end

function internal_force(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  R = create_fields(domain)
  internal_force!(R, U, domain, Uu, state, props, X)
  return R
end

function stiffness(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  assembler = domain.assembler
  stiffness!(assembler, U, domain, Uu, state, props, X)
  K = SparseArrays.sparse!(assembler)
  return 0.5 * (K + K')
end

function internal_force_and_stiffness(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  assembler = domain.assembler
  internal_force_and_stiffness!(assembler, U, domain, Uu, state, props, X)
  # R = @views assembler.residuals[domain.dof.unknown_dofs]
  K = SparseArrays.sparse!(assembler)
  return assembler.residuals, 0.5 * (K + K')
end
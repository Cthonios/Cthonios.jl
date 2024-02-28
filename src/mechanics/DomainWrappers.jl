function energy!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy!(Π, U, sections, state, props, X)

  return nothing
end

function internal_force!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force!(f, U, sections, state, props, X)

  return nothing
end

function stiffness!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  asm = solver.linear_solver.assembler
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  stiffness!(asm, U, sections, state, props, X)

  return nothing
end

function stiffness_action!(Hv::V1, domain::QuasiStaticDomain, Uu, Vv) where V1 <: AbstractArray

  # zero out stuff
  Hv .= zero(eltype(Hv))

  # unpack stuff
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π, V = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  FiniteElementContainers.update_fields!(V, domain.dof, Vv)
  stiffness_action!(Hv, U, sections, state, props, X, V)

  return nothing
end

function energy_and_internal_force!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy_and_internal_force!(Π, f, U, sections, state, props, X)

  return nothing
end

function internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  asm = solver.assembler
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force_and_stiffness!(f, asm, U, sections, state, props, X)

  return nothing
end

function energy_internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  asm = solver.assembler
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy_internal_force_and_stiffness!(Π, f, asm, U, sections, state, props, X)

  return nothing
end

# TODO maybe this one isn't so general...
function energy_internal_force_and_stiffness_action!(solver, domain::QuasiStaticDomain, Uu, Vv)

  # unpack stuff
  Hv = solver.assembler.stiffness_actions
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π, V = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  update_unknowns!(V, domain, Vv)
  energy_internal_force_and_stiffness_action!(Π, f, Hv, U, sections, state, props, X, V)

  return nothing
end

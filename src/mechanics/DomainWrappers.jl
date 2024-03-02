# TODO remove solver from the ones that don't need it...

function energy!(solver, domain::QuasiStaticDomain, Uu, backend)

  # unpack stuff
  sections = domain.sections
  @unpack X, U, state, props, Π, Πs = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy!(Πs, sections, U, state, props, X, backend)
  Π[1] = sum(Πs)

  return nothing
end

function internal_force!(solver, domain::QuasiStaticDomain, Uu, backend)

  # unpack stuff
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force!(f, sections, U, state, props, X, backend)

  return nothing
end

function stiffness!(solver, domain::QuasiStaticDomain, Uu, backend)

  # unpack stuff
  asm = solver.linear_solver.assembler
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  stiffness!(asm, sections, U, state, props, X, backend)

  return nothing
end

# TODO cleanup this method output
function stiffness_action!(Hv::V1, domain::QuasiStaticDomain, Uu, Vv, backend) where V1 <: AbstractArray

  # zero out stuff
  Hv .= zero(eltype(Hv))

  # unpack stuff
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π, V = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  FiniteElementContainers.update_fields!(V, domain.dof, Vv)
  stiffness_action!(Hv, sections, U, state, props, X, V, backend)

  return nothing
end

function energy_and_internal_force!(solver, domain::QuasiStaticDomain, Uu, backend)

  # unpack stuff
  sections = domain.sections
  @unpack X, U, state, props, Π, Πs, f = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy_and_internal_force!(Πs, f, sections, U, state, props, X, backend)
  Π[1] = sum(Πs)
  return nothing
end

function internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, Uu, backend)

  # unpack stuff
  asm = solver.assembler
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force_and_stiffness!(f, asm, sections, U, state, props, X, backend)

  return nothing
end

function energy_internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, Uu, backend)

  # unpack stuff
  asm = solver.assembler
  sections = domain.sections
  @unpack X, U, state, props, Π, Πs, f = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy_internal_force_and_stiffness!(Πs, f, asm, sections, U, state, props, X, backend)
  Π[1] = sum(Πs)
  return nothing
end

# TODO maybe this one isn't so general...
function energy_internal_force_and_stiffness_action!(solver, domain::QuasiStaticDomain, Uu, Vv, backend)

  # unpack stuff
  Hv = solver.assembler.stiffness_actions
  X = domain.domain_cache.X
  f = domain.domain_cache.f
  sections = domain.sections
  @unpack U, state, props, Π, Πs, V = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  update_unknowns!(V, domain, Vv)
  energy_internal_force_and_stiffness_action!(Πs, f, Hv, sections, U, state, props, X, V, backend)
  Π[1] = sum(Πs)
  return nothing
end

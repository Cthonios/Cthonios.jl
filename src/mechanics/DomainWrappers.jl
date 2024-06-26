# TODO remove solver from the ones that don't need it...

function strain_energy!(solver, domain::QuasiStaticDomain, cache, Uu, backend)

  # unpack stuff
  Δt = domain.domain_cache.time.Δt
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, Πs = cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  strain_energy!(Πs, state_new, sections, Δt, X, U, props, state_old, backend)
  Π[1] = sum(Πs)

  return nothing
end

function internal_force!(solver, domain::QuasiStaticDomain, cache, Uu, backend)

  # unpack stuff
  Δt = domain.domain_cache.time.Δt
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, f = cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force!(f, state_new, sections, Δt, X, U, props, state_old, backend)
  
  solver.linear_solver.assembler.residuals .= f

  return nothing
end

function stiffness!(solver, domain::QuasiStaticDomain, cache, Uu, backend)

  # unpack stuff
  asm = solver.linear_solver.assembler
  Δt = domain.domain_cache.time.Δt
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π = cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  stiffness!(asm, state_new, sections, Δt, X, U, props, state_old, backend)

  return nothing
end

# TODO cleanup this method output
function stiffness_action!(solver, domain::QuasiStaticDomain, cache, Uu, Vv, backend)

  # unpack stuff
  Δt = domain.domain_cache.time.Δt
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, V, Hv = cache

  # zero out stuff
  Hv .= zero(eltype(Hv))

  # update fields here
  update_fields!(U, domain, X, Uu)
  FiniteElementContainers.update_fields!(V, domain.dof, Vv)
  stiffness_action!(Hv, state_new, sections, Δt, X, U, props, state_old, V, backend)

  return nothing
end

function strain_energy_and_internal_force!(solver, domain::QuasiStaticDomain, cache, Uu, backend)

  # unpack stuff
  Δt = domain.domain_cache.time.Δt
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, Πs, f = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  strain_energy_and_internal_force!(Πs, f, state_new, sections, Δt, X, U, props, state_old, backend)
  Π[1] = sum(Πs)
  solver.linear_solver.assembler.residuals .= f
  return nothing
end

function internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, cache, Uu, backend)

  # unpack stuff
  Δt = domain.domain_cache.time.Δt
  asm = solver.linear_solver.assembler
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, f = cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force_and_stiffness!(f, asm, state_new, sections, Δt, X, U, props, state_old, backend)
  solver.linear_solver.assembler.residuals .= f

  return nothing
end

function strain_energy_internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, cache, Uu, backend)

  # unpack stuff
  Δt = domain.domain_cache.time.Δt
  asm = solver.linear_solver.assembler
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, Πs, f = cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  strain_energy_internal_force_and_stiffness!(Πs, f, asm, state_new, sections, Δt, X, U, props, state_old, backend)
  Π[1] = sum(Πs)
  solver.linear_solver.assembler.residuals .= f
  return nothing
end

# TODO maybe this one isn't so general...
function strain_energy_internal_force_and_stiffness_action!(solver, domain::QuasiStaticDomain, cache, Uu, Vv, backend)

  # unpack stuff
  Δt = domain.domain_cache.time.Δt
  Hv = solver.linear_solver.assembler.stiffness_actions
  sections = domain.sections
  @unpack X, U, props, state_old, state_new, Π, Πs, f, V = cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  update_unknowns!(V, domain, Vv)
  strain_energy_internal_force_and_stiffness_action!(Πs, f, Hv, state_new, sections, Δt, X, U, props, state_new, V, backend)
  Π[1] = sum(Πs)
  solver.linear_solver.assembler.residuals .= f
  return nothing
end

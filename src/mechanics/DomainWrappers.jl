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
  R = solver.linear_solver.assembler.residuals
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force!(R, U, sections, state, props, X)

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

# function stiffness_action!(solver::AbstractLinearSolver, domain::QuasiStaticDomain, Uu, Vv)

#   # unpack stuff
#   Kv = solver.assembler.stiffness_actions
#   X = domain.domain_cache.X
#   sections = domain.sections
#   @unpack U, state, props, Π, V = domain.domain_cache

#   # update fields here
#   update_fields!(U, domain, X, Uu)
#   update_unknowns!(V, domain, Vv)
#   stiffness_action!(Kv, U, sections, state, props, X, V)

#   return nothing
# end

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
  asm = solver.assembler
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy_and_internal_force!(Π, asm, U, sections, state, props, X)

  return nothing
end

function internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  asm = solver.assembler
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  internal_force_and_stiffness!(asm, U, sections, state, props, X)

  return nothing
end

function energy_internal_force_and_stiffness!(solver, domain::QuasiStaticDomain, Uu)

  # unpack stuff
  asm = solver.assembler
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  energy_internal_force_and_stiffness!(Π, asm, U, sections, state, props, X)

  return nothing
end

function energy_internal_force_and_stiffness_action!(solver, domain::QuasiStaticDomain, Uu, Vv)

  # unpack stuff
  # asm = solver.assembler
  R = solver.assemble.residual
  Hv = solver.assembler.stiffness_actions
  X = domain.domain_cache.X
  sections = domain.sections
  @unpack U, state, props, Π, V = domain.domain_cache

  # update fields here
  update_fields!(U, domain, X, Uu)
  update_unknowns!(V, domain, Vv)
  energy_internal_force_and_stiffness_action!(Π, R, Hv, U, sections, state, props, X, V)

  return nothing
end

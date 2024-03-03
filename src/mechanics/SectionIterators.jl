# top level in place method that loops over sections
# function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X)
function energy!(Πs, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  Πs .= zero(eltype(Πs))
  for (name, section) in pairs(sections)
    Πs_temp        = @views Πs[name]
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    energy!(
      Πs_temp, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, backend
    )
  end 
  return nothing
end 

function internal_force!(f, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    internal_force!(
      f, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, backend
    )
  end 
  return nothing
end 

function stiffness!(assembler::StaticAssembler, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  # update_fields!(U, domain, X, Uu)
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    stiffness!(
      assembler, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, block_count, backend
    )
    block_count = block_count + 1
  end 
  return nothing
end 

function stiffness_action!(Kv, state_new, sections, Δt, X, U, props, state_old, V, backend::Backend)
  Kv .= zero(eltype(Kv))
  for (name, section) in pairs(sections)
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    stiffness_action!(
      Kv, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, V, backend
    )
  end 
  return nothing
end

# dual return outputs
function energy_and_internal_force!(Πs, f, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)
    Πs_temp    = @views Πs[name]
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    energy_and_internal_force!(
      Πs_temp, f, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, backend
    )
  end
  return nothing
end

function internal_force_and_stiffness!(f, assembler, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    internal_force_and_stiffness!(
      f, assembler, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, block_count, backend
    )
    block_count = block_count + 1
  end 
  return nothing
end 


# triple return outputs
function energy_internal_force_and_stiffness!(Πs, f, assembler, state_new, sections, Δt, X, U, props, state_old, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    Πs_temp    = @views Πs[name]
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    energy_internal_force_and_stiffness!(
      Πs_temp, f, assembler, state_new_temp, 
      section, Δt, X, U, props_temp, state_old_temp, block_count, backend
    )
    block_count = block_count + 1
  end
  return nothing
end

function energy_internal_force_and_stiffness_action!(Πs, f, Hv, state_new, sections, Δt, X, U, props, state_old, V, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  Hv .= zero(eltype(Hv))
  for (name, section) in pairs(sections)
    Πs_temp    = @views Πs[name]
    state_old_temp = @views state_old[name]
    state_new_temp = @views state_new[name]
    props_temp     = @views props[name]
    energy_internal_force_and_stiffness_action!(
      Πs_temp, f, Hv, state_new_temp,
      section, Δt, X, U, props_temp, state_old_temp, V, backend
    )
  end
  return nothing
end

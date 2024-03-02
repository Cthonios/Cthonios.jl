# top level in place method that loops over sections
# function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X)
function energy!(Πs, sections, U, state, props, X, backend::Backend)
  Πs .= zero(eltype(Πs))
  for (name, section) in pairs(sections)
    Πs_temp    = @views Πs[name]
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy!(Πs_temp, section, U, state_temp, props_temp, X, backend)
  end 
  return nothing
end 

function internal_force!(f, sections, U, state, props, X, backend::Backend)
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    internal_force!(f, section, U, state_temp, props_temp, X, backend)
  end 
  return nothing
end 

function stiffness!(assembler::StaticAssembler, sections, U, state, props, X, backend::Backend)
  # update_fields!(U, domain, X, Uu)
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    stiffness!(assembler, section, U, state_temp, props_temp, X, block_count, backend)
    block_count = block_count + 1
  end 
  return nothing
end 

function stiffness_action!(Kv, sections, U, state, props, X, V, backend::Backend)
  Kv .= zero(eltype(Kv))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    stiffness_action!(Kv, section, U, state_temp, props_temp, X, V, backend)
  end 
  return nothing
end

# dual return outputs
function energy_and_internal_force!(Πs, f, sections, U, state, props, X, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)
    Πs_temp    = @views Πs[name]
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_and_internal_force!(Πs_temp, f, section, U, state_temp, props_temp, X, backend)
  end
  return nothing
end

function internal_force_and_stiffness!(f, assembler, sections, U, state, props, X, backend::Backend)
  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    internal_force_and_stiffness!(f, assembler, section, U, state_temp, props_temp, X, block_count, backend)
    block_count = block_count + 1
  end 
  return nothing
end 


# triple return outputs
function energy_internal_force_and_stiffness!(Πs, f, assembler, sections, U, state, props, X, backend::Backend)
  Πs .= zero(eltype(Πs))
  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    Πs_temp    = @views Πs[name]
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_internal_force_and_stiffness!(Πs_temp, f, assembler, section, U, state_temp, props_temp, X, block_count, backend)
    block_count = block_count + 1
  end
  return nothing
end

function energy_internal_force_and_stiffness_action!(Πs, f, Hv, sections, U, state, props, X, V, backend::Backend)
  Πs .= zero(eltype(Π))
  f .= zero(eltype(f))
  Hv .= zero(eltype(Hv))
  for (name, section) in pairs(sections)
    Πs_temp    = @views Πs[name]
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_internal_force_and_stiffness_action!(Πs_temp, f, Hv, section, U, state_temp, props_temp, X, V, backend)
  end
  return nothing
end

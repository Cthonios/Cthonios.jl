# top level in place method that loops over sections
# function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X)
function energy!(Π, U, sections, state, props, X)
  Π .= zero(eltype(Π))
  # update_fields!(U, domain, X, Uu)
  # for (name, section) in pairs(domain.sections)
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    Π[1] = Π[1] + energy(section, U, state_temp, props_temp, X)
    # Π .+= energy(section, U, state_temp, props_temp, X)
  end 
  return nothing
end 

function internal_force!(R, U, sections, state, props, X)
  R .= zero(eltype(R))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    internal_force!(R, section, U, state_temp, props_temp, X)
  end 
  return nothing
end 

function stiffness!(assembler::StaticAssembler, U, sections, state, props, X)
  # update_fields!(U, domain, X, Uu)
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    stiffness!(assembler, section, U, state_temp, props_temp, X, block_count)
    block_count = block_count + 1
  end 
  return nothing
end 

function stiffness_action!(Kv, U, sections, state, props, X, V)
  Kv .= zero(eltype(Kv))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    stiffness_action!(Kv, section, U, state_temp, props_temp, X, V)
  end 
  return nothing
end

# multi return outputs
function internal_force_and_stiffness!(assembler, U, sections, state, props, X)
  # kernel setup, TODO add backend as argument
  backend = CPU()
  kernel! = internal_force_and_stiffness_kernel!(backend)

  assembler.residuals .= zero(eltype(assembler.residuals))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    kernel!(assembler, section, U, state_temp, props_temp, X, block_count, ndrange=(FiniteElementContainers.num_q_points(section), num_elements(section)))
    block_count = block_count + 1
  end 
  return nothing
end 

function energy_and_internal_force!(Π, assembler, U, sections, state, props, X)
  Π .= zero(eltype(Π))
  assembler.residuals .= zero(eltype(assembler.residuals))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_and_internal_force!(Π, assembler, section, U, state_temp, props_temp, X)
  end
  return nothing
end

function energy_internal_force_and_stiffness!(Π, assembler, U, sections, state, props, X)
  Π .= zero(eltype(Π))
  assembler.residuals .= zero(eltype(assembler.residuals))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_internal_force_and_stiffness!(Π, assembler, section, U, state_temp, props_temp, X, block_count)
    block_count = block_count + 1
  end
  return nothing
end

function energy_internal_force_and_stiffness_action!(Π, R, Hv, U, sections, state, props, X, V)
  Π .= zero(eltype(Π))
  R .= zero(eltype(R))
  Hv .= zero(eltype(Hv))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_internal_force_and_stiffness_action!(Π, R, Hv, section, U, state_temp, props_temp, X, V)
  end
  return nothing
end

# top level in place method that loops over sections
# function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X)
function energy!(Π, sections, U, state, props, X, backend::NoBackend)
  Π .= zero(eltype(Π))
  # update_fields!(U, domain, X, Uu)
  # for (name, section) in pairs(domain.sections)
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    Π[1] = Π[1] + energy(section, U, state_temp, props_temp, X, backend)
    # Π .+= energy(section, U, state_temp, props_temp, X)
  end 
  return nothing
end 

function energy!(Π, Πs, sections, U, state, props, X, backend::CPU)

  # kernel setup
  kernel! = energy_kernel!(backend)

  for (name, section) in pairs(sections)
    NQ = FiniteElementContainers.num_q_points(section)
    NE = num_elements(section)
    Πs_temp    = @views Πs[name]
    state_temp = @views state[name]
    props_temp = @views props[name]
    kernel!(
      Πs_temp, 
      section,
      U, state_temp, props_temp, X,
      ndrange=(NQ, NE)
    )
  end 
  Π[1] = sum(Πs)
  # Π[1] = 0.0
  return nothing
end 

function kernel_iterator!(
  assembler::StaticAssembler,  # outputs
  kernel!,                     # kernel to dispatch on
  sections, U, state, props, X # inputs to kernel
)
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    NQ = FiniteElementContainers.num_q_points(section)
    NE = num_elements(section)
    state_temp = @views state[name]
    props_temp = @views props[name]
    kernel!(
      assembler, # outputs of kernel 
      section, U, state_temp, props_temp, X, block_count, # inputs of kernel
      ndrange=(NQ, NE)
    )
    block_count = block_count + 1
  end
end

function internal_force!(f, sections, U, state, props, X, backend::NoBackend)
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    internal_force!(f, section, U, state_temp, props_temp, X, backend)
  end 
  return nothing
end 

function stiffness!(assembler::StaticAssembler, sections, U, state, props, X, backend::NoBackend)
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

function stiffness_action!(Kv, sections, U, state, props, X, V, backend::NoBackend)
  Kv .= zero(eltype(Kv))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    stiffness_action!(Kv, section, U, state_temp, props_temp, X, V, backend)
  end 
  return nothing
end

# multi return outputs
function internal_force_and_stiffness!(f, assembler, sections, U, state, props, X, backend::NoBackend)
  # kernel setup, TODO add backend as argument
  # backend = CPU()
  # kernel! = internal_force_and_stiffness_kernel!(backend)

  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    # kernel!(f, assembler, section, U, state_temp, props_temp, X, block_count, ndrange=(FiniteElementContainers.num_q_points(section), num_elements(section)))
    internal_force_and_stiffness!(f, assembler, section, U, state_temp, props_temp, X, block_count, backend)
    block_count = block_count + 1
  end 
  return nothing
end 

function energy_and_internal_force!(Π, f, sections, U, state, props, X, backend::NoBackend)
  Π .= zero(eltype(Π))
  f .= zero(eltype(f))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_and_internal_force!(Π, f, section, U, state_temp, props_temp, X, backend)
  end
  return nothing
end

function energy_internal_force_and_stiffness!(Π, f, assembler, sections, U, state, props, X, backend::NoBackend)
  Π .= zero(eltype(Π))
  f .= zero(eltype(f))
  assembler.stiffnesses .= zero(eltype(assembler.stiffnesses))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_internal_force_and_stiffness!(Π, f, assembler, section, U, state_temp, props_temp, X, block_count, backend)
    block_count = block_count + 1
  end
  return nothing
end

function energy_internal_force_and_stiffness_action!(Π, f, Hv, sections, U, state, props, X, V, backend::NoBackend)
  Π .= zero(eltype(Π))
  f .= zero(eltype(f))
  Hv .= zero(eltype(Hv))
  for (name, section) in pairs(sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    energy_internal_force_and_stiffness_action!(Π, f, Hv, section, U, state_temp, props_temp, X, V, backend)
  end
  return nothing
end

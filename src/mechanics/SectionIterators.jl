# top level in place method that loops over sections
function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X)
  update_fields!(U, domain, X, Uu)
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    Π[1] = Π[1] + energy(section, U, state_temp, props_temp, X)
    # Π .+= energy(section, U, state_temp, props_temp, X)
  end 
  return nothing
end 

function internal_force!(R, U, domain::QuasiStaticDomain, Uu, state, props, X)
  update_fields!(U, domain, X, Uu)
  R .= zero(eltype(R))
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    internal_force!(R, section, U, state_temp, props_temp, X)
  end 
  return nothing
end 

function stiffness!(assembler, U, domain::QuasiStaticDomain, Uu, state, props, X)
  update_fields!(U, domain, X, Uu)
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    stiffness!(assembler, section, U, state_temp, props_temp, X, block_count)
    block_count = block_count + 1
  end 
  return nothing
end 

# multi return outputs
function internal_force_and_stiffness!(assembler, U, domain::QuasiStaticDomain, Uu, state, props, X)
  backend = CPU()
  kernel! = internal_force_and_stiffness_kernel!(backend)

  update_fields!(U, domain, X, Uu)
  assembler.residuals .= zero(eltype(assembler.residuals))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    kernel!(assembler, section, U, state_temp, props_temp, X, block_count, ndrange=(FiniteElementContainers.num_q_points(section), num_elements(section)))
    block_count = block_count + 1
  end 
  return nothing
end 

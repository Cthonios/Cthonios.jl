# TODO have the volume integral own
# physics methods for energy, gradient, hessian

struct VolumeIntegral{F1, F2, F3, F4} <: AbstractIntegral
  integrand::F1
  integrand_gradient_u::F2
  integrand_gradient_x::F3
  integrand_hessian_u::F4
end

# constructor for energy based physics
function VolumeIntegral(::typeof(energy))
  # TODO change to have use AD as a flag to flag whether or not to race energy
  # wrt certain things like u that could easily be defined analytically
  # grad_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> energy(phys, cell, z, x, state, props, t), u)
  grad_func = gradient
  grad_x_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> energy(phys, cell, u, z, state, props, t), x)
  # hess_func(phys, cell, u, x, state, props, t) = ForwardDiff.hessian(z -> energy(phys, cell, z, x, state, props, t), u)
  hess_func = hessian
  return VolumeIntegral(energy, grad_func, grad_x_func, hess_func)
end

"""
$(TYPEDSIGNATURES)
Volume integral over a over a domain ```domain``` to fill a global value
```global_val``` based on a quadrature level function 
```f``` provided a nodal field ```U``` and set of
paramaters ```p```. This method is useful for filling
quantities such as objectives, gradients, or hessians.
"""
function integrate!(global_val, ::VolumeIntegral, func, Uu, p, domain)
  update_field_bcs!(current_solution(p), domain, p.Ubc)
  update_field_unknowns!(current_solution(p), domain, Uu)

  # loop over sections
  for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
    ND, NN, NP, NS = size(section)
    fspace, physics = section.fspace, section.physics
    # loop over elements
    for e in 1:num_elements(fspace)
      X_el = element_coordinates(section, p.X, e)
      U_el = element_fields(section, current_solution(p), e)
      # props_el = element_props(section, p.props, section_name, e)
      props_el = @views SVector{NP, eltype(p.props)}(p.props[section_name])
      local_val = scratch_variable(global_val, section)
      # loop over quadrature points
      for q in 1:num_q_points(fspace)
        # interps = MappedInterpolants(fspace.ref_fe, X_el, q)
        interps = fspace.ref_fe.cell_interps.vals[q]
        state = @views SVector{num_states(section.physics), Float64}(p.state_old[section_name][:, q, e])
        local_val += func(physics, interps, U_el, X_el, state, props_el, p.t)
      end

      # assembly
      assemble!(global_val, fspace, block_num, e, local_val)
    end
  end

  return nothing
end 

"""
$(TYPEDSIGNATURES)
Volume integral over a domain ```domain``` to fill a global value
```global_val``` based on a quadrature level function 
```f``` provided a nodal field ```U```, a set of
paramaters ```p```, and a vector ```V```. This method
is useful for quantities such as hessian vector productions.
"""
function integrate!(global_val, ::VolumeIntegral, func, Uu, p, Vv, domain)
  # update_field_bcs!(current_solution(p), domain, current_solution(p)bc)
  # update_field_unknowns!(current_solution(p), domain, Uu)
  # update_field_unknowns!(p.hvp_scratch, domain, Vv)
  update_field_bcs!(current_solution(p), domain, p.Ubc)
  update_field_unknowns!(current_solution(p), domain, Uu)
  update_field_unknowns!(p.hvp_scratch, domain, Vv)

  # update_field_bcs!(U, domain, Ubc)
  for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
    ND, NN, NP, NS = size(section)
    fspace, physics = section.fspace, section.physics
    
    # loop over elements
    for e in 1:num_elements(fspace)
      ND, NN, NP, NS = size(section)
      NF = num_dofs_per_node(section.fspace)
      X_el = element_coordinates(section, p.X, e)
      U_el = element_fields(section, current_solution(p), e)
      V_el = element_fields(section, p.hvp_scratch, e)
      props_el = @views SVector{NP, eltype(p.props)}(p.props[section_name])
      local_val = zeros(SMatrix{NN * NF, NN * NF, eltype(global_val), NN * NF * NN * NF})
      # loop over quadrature points
      for q in 1:num_q_points(fspace)
        # TODO need to modify this getindex. We can't trace this
        # interps = MappedInterpolants(fspace.ref_fe, X_el, q)
        interps = fspace.ref_fe.cell_interps.vals[q]
        state = @views SVector{num_states(section.physics), Float64}(p.state_old[section_name][:, q, e])
        local_val += func(physics, interps, U_el, X_el, state, props_el, p.t)
      end
      local_val = local_val * vec(V_el)

      # assembly
      assemble!(global_val, fspace, block_num, e, local_val)
    end
  end

  # shouldn't need any neumann stuff here

  return nothing
end

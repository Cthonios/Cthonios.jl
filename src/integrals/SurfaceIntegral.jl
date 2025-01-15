struct SurfaceIntegral{F1, F2, F3, F4} <: AbstractIntegral
  integrand::F1
  integrand_gradient_u::F2
  integrand_gradient_x::F3
  integrand_hessian_u::F4
end

function SurfaceIntegral(::typeof(neumann_energy))
  # TODO change to have use AD as a flag to flag whether or not to race energy
  # wrt certain things like u that could easily be defined analytically
  # grad_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> energy(phys, cell, z, x, state, props, t), u)
  grad_func = neumann_gradient
  # grad_x_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> neumann_energy(phys, cell, u, z, state, props, t), x)
  grad_x_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_energy(phys, cell, u, z, t, bc, val, r, q, f), x)

  # hess_func(phys, cell, u, x, state, props, t) = ForwardDiff.hessian(z -> energy(phys, cell, z, x, state, props, t), u)
  hess_func = neumann_hessian
  return SurfaceIntegral(energy, grad_func, grad_x_func, hess_func)
end

function integrate!(global_val, ::SurfaceIntegral, func, Uu, p, domain::Domain)
  update_field_bcs!(current_solution(p), domain, p.Ubc)
  update_field_unknowns!(current_solution(p), domain, Uu)

  # Neumann BCs
  # TODO minor inefficiency for hessian by 
  # assembling a zero matrix for each side a traction
  # acts on. Maybe the compiler is smart enough to compile this away?
  # Maybe we should break this method up?
  bc_index = 1
  for (block_num, (section_name, section)) in enumerate(pairs(domain.neumann_bc_sections))
    ND, NN, NP, NS = size(section)
    bc, fspace, physics = section.bc, section.fspace, section.physics
    for e in 1:num_elements(fspace)
      X_el = surface_element_coordinates(section, p.X, e)
      U_el = surface_element_fields(section, current_solution(p), e)
      local_val = surface_scratch_variable(global_val, section)
      for q in 1:ReferenceFiniteElements.num_quadrature_points(fspace.ref_fe.surface_element)
        bc_val = p.nbc[bc_index]
        face = section.bc.sides[e]
        interps = fspace.ref_fe.surface_interps.vals[q, face]
        local_val += func(physics, interps, U_el, X_el, p.t, bc, bc_val, fspace.ref_fe, q, face)
        bc_index = bc_index + 1
      end

      # assembly
      assemble!(global_val, fspace, block_num, e, local_val)
    end
  end
  return nothing
end

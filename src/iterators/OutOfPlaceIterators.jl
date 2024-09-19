function domain_iterator(f, domain, Uu, p)
  # setup field
  U = create_fields(domain)
  U[domain.dirichlet_dofs] = p.Ubc
  U[domain.dof.unknown_dofs] = Uu

  # loop over all sections/elements/q points
  value = zero(eltype(U))
  for sec in domain.sections
    for e in 1:num_elements(sec.fspace)
      dof_conn = dof_connectivity(sec.fspace, e)
      U_el = element_fields(sec, U, dof_conn)
      for q in 1:num_q_points(sec.fspace)
        interps = getindex(sec.fspace, p.X, q, e)
        value += f(sec.physics, interps, U_el)
      end
    end
  end
  return value
end

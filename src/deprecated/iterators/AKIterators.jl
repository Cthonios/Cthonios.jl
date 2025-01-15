# move to FEMContainers eventually
# function AK.Ker
function KA.get_backend(arr::T) where T <: FiniteElementContainers.AbstractField
  return KA.get_backend(arr.vals)
end

# update bcs/unknowns methods prototypes
function update_field_bcs!(U, domain, Ubc, ::AKIteator)
  @assert length(domain.dirichlet_dofs) == length(Ubc)
  AK.foreachindex(domain.dirichlet_dofs) do n
    dof = domain.dirichlet_dofs[n]
    U[dof] = Ubc[n]
  end
  return nothing
end

function update_field_unknowns!(U, domain, Uu, ::AKIteator)
  @assert length(Uu) == length(domain.dof.unknown_indices)
  AK.foreachindex(domain.dof.unknown_indices) do n
    dof = domain.dof.unknown_indices[n]
    U[dof] = Uu[n]
  end
  return nothing
end

function domain_scalar_iterator!(
  global_val, 
  f, 
  domain::Domain, Uu, p::ObjectiveParameters,
  iterator::AKIteator
)
  # update whole field U
  update_field_bcs!(p.U, domain, p.Ubc, iterator)
  update_field_unknowns!(p.U, domain, Uu)

  # loop over sections
  for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
    ND, NN, NP, NS = size(section)
    fspace, physics = section.fspace, section.physics
    temp_vals = p.q_vals_scratch[section_name]
    # for e in 1:num_elements(fspace)
    AK.foraxes(fspace.conn, 2) do e
      X_el = element_coordinates(section, p.X, e)
      U_el = element_fields(section, p.U, e)
      # props_el = element_props(section, p.props, section_name, e)
      props_el = @views SVector{NP, eltype(p.props)}(p.props[section_name])
      # props_el = SVector{NP, eltype(p.props)}((1., 1.))
      # props_el = SVector{2, Float64}((1., 1.))
      local_val = scratch_variable(global_val, section)
      # loop over quadrature points
      for q in 1:num_q_points(fspace)
        # interps = MappedInterpolants(fspace.ref_fe, X_el, q)
        interps = fspace.ref_fe.cell_interps.vals[q]
        # state = @views SVector{num_states(section.physics), Float64}(p.state_old[section_name][:, q, e])
        state = SVector{0, Float64}()
        # Atomix.@atomic local_val[:] += f(physics, interps, U_el, X_el, state, props_el, p.t)
        # p.q_vals_scratch[section_name][q, e] = f(physics, interps, U_el, X_el, state, props_el, p.t)
        temp_vals[q, e] = f(physics, interps, U_el, X_el, state, props_el, p.t)
      end

      # assembly
      # assemble!(global_val, fspace, block_num, e, local_val)
    end
  end
end
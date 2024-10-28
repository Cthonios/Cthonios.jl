"""
$(TYPEDSIGNATURES)
"""
function element_coordinates(section, X, e)
  ND, NN, _, _ = size(section)
  conn = connectivity(section.fspace, e)
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[:, conn])
  return X_el
end

# Try and combine below and above into single method
function surface_element_coordinates(section, X, e)
  ND, _, _, _ = size(section)
  NN = num_nodes_per_side(section)
  conn = connectivity(section.fspace, e)
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[:, conn])
  return X_el
end

"""
$(TYPEDSIGNATURES)
"""
function element_fields(section, U, e)
  _, NN, _, _ = size(section)
  NF = num_dofs_per_node(section.fspace)
  dof_conn = dof_connectivity(section.fspace, e)
  # U_el = SMatrix{NF, NN, eltype(U), NF * NN}(@views U[dof_conn])
  U_el = SVector{NF * NN, eltype(U)}(@views U[dof_conn])
  return U_el
end

function surface_element_fields(section, U, e)
  # _, NN, _, _ = size(section)
  NN = num_nodes_per_side(section)
  NF = num_dofs_per_node(section.fspace)
  dof_conn = dof_connectivity(section.fspace, e)
  # U_el = SMatrix{NF, NN, eltype(U), NF * NN}(@views U[dof_conn])
  U_el = SVector{NF * NN, eltype(U)}(@views U[dof_conn])
  return U_el
end

"""
$(TYPEDSIGNATURES)
Setup a scratch variable for an energy
like calculation
"""
function scratch_variable(global_val::Vector, section)
  return zero(eltype(global_val))
end

function surface_scratch_variable(global_val::Vector, section)
  return zero(eltype(global_val))
end

"""
$(TYPEDSIGNATURES)
Setup a scratch variable for a force
like calculation
"""
function scratch_variable(global_val::T, section) where T <: Union{Matrix, NodalField}
  ND, NN, NP, NS = size(section)
  NF = num_dofs_per_node(section.fspace)
  # TODO change to SVector
  # return zeros(SMatrix{NF, NN, eltype(global_val), NF * NN})
  return zeros(SVector{NF * NN, eltype(global_val)})
end

function surface_scratch_variable(global_val::T, section) where T <: Union{Matrix, NodalField}
  # ND, NN, NP, NS = size(section)
  NN = num_nodes_per_side(section)
  NF = num_dofs_per_node(section.fspace)
  # TODO change to SVector
  # return zeros(SMatrix{NF, NN, eltype(global_val), NF * NN})
  return zeros(SVector{NF * NN, eltype(global_val)})
end

"""
$(TYPEDSIGNATURES)
Setup a scratch variable for a stiffness
like calculation
"""
function scratch_variable(global_val::FiniteElementContainers.StaticAssembler, section)
  ND, NN, NP, NS = size(section)
  NF = num_dofs_per_node(section.fspace)
  return zeros(SMatrix{NF * NN, NF * NN, eltype(global_val.stiffnesses), NF * NN * NF * NN})
end

function surface_scratch_variable(global_val::FiniteElementContainers.StaticAssembler, section)
  # ND, NN, NP, NS = size(section)
  NN = num_nodes_per_side(section)
  NF = num_dofs_per_node(section.fspace)
  return zeros(SMatrix{NF * NN, NF * NN, eltype(global_val.stiffnesses), NF * NN * NF * NN})
end

"""
$(TYPEDSIGNATURES)
Iterator over a domain ```domain``` to fill a global value
```global_val``` based on a quadrature level function 
```f``` provided a nodal field ```U``` and set of
paramaters ```p```. This method is useful for filling
quantities such as objectives, gradients, or hessians.
"""
function domain_iterator!(global_val, f, domain::Domain, Uu, p::ObjectiveParameters)
  update_field_bcs!(p.U, domain, p.Ubc)
  update_field_unknowns!(p.U, domain, Uu)

  # loop over sections
  for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
    ND, NN, NP, NS = size(section)
    fspace, physics = section.fspace, section.physics
    # loop over elements
    for e in 1:num_elements(fspace)
      X_el = element_coordinates(section, p.X, e)
      U_el = element_fields(section, p.U, e)
      # props_el = element_props(section, p.props, section_name, e)
      props_el = @views SVector{NP, eltype(p.props)}(p.props[section_name])
      local_val = scratch_variable(global_val, section)
      # loop over quadrature points
      for q in 1:num_q_points(fspace)
        interps = MappedInterpolants(fspace.ref_fe, X_el, q)
        state = @views SVector{num_states(section.physics), Float64}(p.state_old[section_name][:, q, e]...)
        local_val += f(physics, interps, U_el, p.t, state, props_el)
      end

      # assembly
      assemble!(global_val, fspace, block_num, e, local_val)
    end
  end

  # Neumann BCs
  # TODO minor inefficiency for hessian by 
  # assembling a zero matrix for each side a traction
  # acts on. Maybe the compiler is smart enough to compile this away?
  # Maybe we should break this method up?
  for (block_num, (section_name, section)) in enumerate(pairs(domain.neumann_bc_sections))
    ND, NN, NP, NS = size(section)
    bc, fspace, physics = section.bc, section.fspace, section.physics
    for e in 1:num_elements(fspace)
      X_el = surface_element_coordinates(section, p.X, e)
      U_el = surface_element_fields(section, p.U, e)
      # props_el = nothing
      local_val = scratch_variable(global_val, section)
      # for q in 1:num_q_points(fspace)
      for q in 1:ReferenceFiniteElements.num_quadrature_points(fspace.ref_fe.surface_element)
        face = section.bc.sides[e]
        interps = MappedSurfaceInterpolants(fspace.ref_fe, X_el, q, face)
        # TODO what to do about state? Do we interpolate or just never pass it?
        # state = nothing
        # @show size(local_val)
        local_val += f(physics, bc, interps, U_el, p.t)
      end

      # assembly
      assemble!(global_val, fspace, block_num, e, local_val)
    end
  end
  return nothing
end

"""
$(TYPEDSIGNATURES)
Iterator over a domain ```domain``` to fill a global value
```global_val``` based on a quadrature level function 
```f``` provided a nodal field ```U```, a set of
paramaters ```p```, and a vector ```V```. This method
is useful for quantities such as hessian vector productions.
"""
# function domain_iterator!(global_val, f, domain::Domain, U::NodalField, X::NodalField, V::NodalField)
function domain_iterator!(global_val, f, domain::Domain, Uu, p::ObjectiveParameters, Vv)
  update_field_bcs!(p.U, domain, p.Ubc)
  update_field_unknowns!(p.U, domain, Uu)
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
      U_el = element_fields(section, p.U, e)
      V_el = element_fields(section, p.hvp_scratch, e)
      props_el = @views SVector{NP, eltype(p.props)}(p.props[section_name])
      local_val = zeros(SMatrix{NN * NF, NN * NF, eltype(global_val), NN * NF * NN * NF})
      # loop over quadrature points
      for q in 1:num_q_points(fspace)
        # TODO need to modify this getindex. We can't trace this
        interps = MappedInterpolants(fspace.ref_fe, X_el, q)
        state = @views SVector{num_states(section.physics), Float64}(p.state_old[section_name][:, q, e]...)
        local_val += f(physics, interps, U_el, p.t, state, props_el)
      end
      local_val = local_val * vec(V_el)

      # assembly
      assemble!(global_val, fspace, block_num, e, local_val)
    end
  end

  # shouldn't need any neumann stuff here

  return nothing
end

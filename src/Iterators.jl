"""
$(TYPEDSIGNATURES)
"""
function element_coordinates(section, X, conn)
  ND, NN, _, _ = size(section)
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[:, conn])
  return X_el
end

"""
$(TYPEDSIGNATURES)
"""
function element_fields(section, U, dof_conn)
  _, NN, _, _ = size(section)
  # NF = FiniteElementContainers.num_fields(U)
  NF = num_dofs_per_node(section.fspace)
  U_el = SMatrix{NF, NN, eltype(U), NF * NN}(@views U[dof_conn])
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

"""
$(TYPEDSIGNATURES)
Iterator over a domain ```domain``` to fill a global value
```global_val``` based on a quadrature level function 
```f``` provided a nodal field ```U``` and set of
paramaters ```p```. This method is useful for filling
quantities such as objectives, gradients, or hessians.
"""
function domain_iterator!(global_val, f, domain, Uu, p::ObjectiveParameters)
  update_field_unknowns!(p.U, domain, Uu)
  domain_iterator!(global_val, f, domain, p.U, p.Ubc, p.X)
  return nothing
end

function domain_iterator!(global_val, f, domain::Domain, U::NodalField, Ubc, X::NodalField)
  update_field_bcs!(U, domain, Ubc)
  domain_iterator!(global_val, f, domain, U, X)
  return nothing
end

"""
$(TYPEDSIGNATURES)
Iterator over a domain ```domain``` to fill a global value
```global_val``` based on a quadrature level function 
```f``` provided a nodal field ```U``` and set of
paramaters ```p```. This method is useful for filling
quantities such as objectives, gradients, or hessians.
"""
function domain_iterator!(global_val, f, domain::Domain, U, X::NodalField)
  # loop over sections
  for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
    fspace, physics = section.fspace, section.physics
    # loop over elements
    for e in 1:num_elements(fspace)
      dof_conn = dof_connectivity(fspace, e)
      U_el = element_fields(section, U, dof_conn)
      local_val = scratch_variable(global_val, section)
      # loop over quadrature points
      for q in 1:num_q_points(fspace)
        interps = getindex(fspace, X, q, e)
        local_val += f(physics, interps, U_el)
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
function domain_iterator!(global_val, f, domain::Domain, Uu, p::ObjectiveParameters, Vv)
  update_field_unknowns!(p.U, domain, Uu)
  update_field_unknowns!(p.hvp_scratch, domain, Vv)
  domain_iterator!(global_val, f, domain, p.U, p.Ubc, p.X, p.hvp_scratch)
  return nothing
end

function domain_iterator!(global_val, f, domain::Domain, U::NodalField, Ubc, X::NodalField, V::NodalField)
  update_field_bcs!(U, domain, Ubc)
  domain_iterator!(global_val, f, domain, U, X, V)
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
function domain_iterator!(global_val, f, domain::Domain, U::NodalField, X::NodalField, V::NodalField)
  # update_field_bcs!(U, domain, Ubc)
  for (block_num, (section_name, section)) in enumerate(pairs(domain.sections))
    fspace, physics = section.fspace, section.physics
    # loop over elements
    for e in 1:num_elements(fspace)
      ND, NN, NP, NS = size(section)
      NF = num_dofs_per_node(section.fspace)
      dof_conn = dof_connectivity(fspace, e)
      U_el = element_fields(section, U, dof_conn)
      V_el = element_fields(section, V, dof_conn)
      local_val = zeros(SMatrix{NN * NF, NN * NF, eltype(global_val), NN * NF * NN * NF})
      # loop over quadrature points
      for q in 1:num_q_points(fspace)
        interps = getindex(fspace, X, q, e)
        local_val += f(physics, interps, U_el)
      end
      local_val = local_val * vec(V_el)

      # assembly
      assemble!(global_val, fspace, block_num, e, local_val)
    end
  end
  return nothing
end

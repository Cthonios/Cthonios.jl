# TODO move some of these helper methods to FiniteElementContainers
"""
$(TYPEDSIGNATURES)
"""
function element_coordinates(section, X, e)
  ND, NN, _, _ = size(section)
  conn = connectivity(section.fspace, e)
  # X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[:, conn])
  X_el = SVector{ND * NN, eltype(X)}(@views X[:, conn])
  return X_el
end

# Try and combine below and above into single method
function surface_element_coordinates(section, X, e)
  ND, _, _, _ = size(section)
  NN = num_nodes_per_side(section)
  conn = connectivity(section.fspace, e)
  # X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[:, conn])
  X_el = SVector{ND * NN, eltype(X)}(@views X[:, conn])
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
  # NN = ReferenceFiniteElements.num_vertices_per_edge(section.fspace.ref_fe)
  NN = ReferenceFiniteElements.num_vertices(section.fspace.ref_fe.surface_element)
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
function scratch_variable(global_val::FiniteElementContainers.Assembler, section)
  ND, NN, NP, NS = size(section)
  NF = num_dofs_per_node(section.fspace)
  return zeros(SMatrix{NF * NN, NF * NN, eltype(global_val.stiffnesses), NF * NN * NF * NN})
end

function surface_scratch_variable(global_val::FiniteElementContainers.Assembler, section)
  # ND, NN, NP, NS = size(section)
  # NN = ReferenceFiniteElements.num_nodes_per_side(section.fspace.ref_fe)
  NN = ReferenceFiniteElements.num_vertices(section.fspace.ref_fe.surface_element)
  NF = num_dofs_per_node(section.fspace)
  return zeros(SMatrix{NF * NN, NF * NN, eltype(global_val.stiffnesses), NF * NN * NF * NN})
end

function unpack_state(section, state, q, e)
  NS = ConstitutiveModels.num_state_vars(section.model)
  state_q = SVector{NS, eltype(state)}(@views state[:, q, e])
  return state_q
end

function unpack_element(section, U, props, X, e)
  # get sizes
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NP = ConstitutiveModels.num_properties(section.model)

  # setup arrays
  conn = dof_connectivity(section, e)
  U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
  props_el = SVector{NP, eltype(props)}(@views props[:, e])
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
  return conn, U_el, props_el, X_el
end

function unpack_element(section, U, props, X, V, e)
  # get sizes
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NP = ConstitutiveModels.num_properties(section.model)
  NDOF = ND * NN

  # setup arrays
  conn = dof_connectivity(section, e)
  U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
  props_el = SVector{NP, eltype(props)}(@views props[:, e])
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
  V_el = SVector{NDOF, eltype(U)}(@views V[conn])
  return conn, U_el, props_el, X_el, V_el
end

function unpack(section, U, state, props, X, q, e)
  # get sizes
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)

  # setup arrays
  conn     = dof_connectivity(section, e)
  U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
  props_el = SVector{NP, eltype(props)}(@views props[:, e])
  state_q = SVector{NS, eltype(state)}(@views state[:, q, e])
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
  return conn, U_el, state_q, props_el, X_el
end

function interpolants(section, q)
  N    = shape_function_values(section, q)
  ∇N_ξ = shape_function_gradients(section, q)
  w    = quadrature_weights(section, q)
  return N, ∇N_ξ, w
end

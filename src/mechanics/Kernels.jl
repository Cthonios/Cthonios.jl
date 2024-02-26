# hot loop over quadrature points for a single section
function energy(section::TotalLagrangeSection, U, state, props, X)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)

  W = 0.0

  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)

    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    props_el = SVector{NP, eltype(props)}(@views props[:, e])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])

      W = W + energy(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
    end
  end

  return W
end

function internal_force!(R, section::TotalLagrangeSection, U, state, props, X)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)

  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)

    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    props_el = SVector{NP, eltype(props)}(@views props[:, e])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    R_el = zeros(SVector{ND * NN, eltype(R)})
    
    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])
      R_q = internal_force(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
      # @show size(temp)
      # @assert false
      R_el = R_el + R_q
    end

    # TODO make a general method in FiniteElementContainers.jl
    for n in 1:ND * NN
      R[conn[n]] += R_el[n]
    end
  end
end

function stiffness!(assembler, section::TotalLagrangeSection, U, state, props, X, block_id)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)

    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    props_el = SVector{NP, eltype(props)}(@views props[:, e])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    # K_el = zeros(SVector{ND * NN, eltype(R)})
    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])
      K_el = K_el + stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
    end

    # assemble into global matrix
    assemble!(assembler, K_el, block_id, e)
  end
end

function stiffness_action!(Kv, section::TotalLagrangeSection, U, state, props, X, V)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)

    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    props_el = SVector{NP, eltype(props)}(@views props[:, e])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    # V_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views V[conn])
    V_el = SVector{NDOF, eltype(U)}(@views V[conn])
    # K_el = zeros(SVector{ND * NN, eltype(R)})
    # K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})
    Kv_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])
      K_q = stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
      Kv_el = Kv_el + K_q * V_el
    end

    # assemble into global matrix
    # assemble!(assembler, K_el, block_id, e)
    assemble!(Kv, Kv_el, conn)
  end
end

@kernel function internal_force_and_stiffness_kernel!(assembler, section::TotalLagrangeSection, U, state, props, X, block_id)
  q, e = @index(Global, NTuple)

  # stuff for static arrays
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)

  # unpack arrays
  conn = dof_connectivity(section, e)

  U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
  props_el = SVector{NP, eltype(props)}(@views props[:, e])
  X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
  state_q = SVector{NS, eltype(props)}(@views state[:, q, e])

  # run routine
  R_q, K_q = internal_force_and_stiffness(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    shape_function_values(section, q),
    shape_function_gradients(section, q),
    quadrature_weights(section, q)
  )

  # update arrays
  FiniteElementContainers.assemble_atomic!(assembler, R_q, conn)
  FiniteElementContainers.assemble_atomic!(assembler, K_q, block_id, e)

  return nothing
end

function energy_and_internal_force!(Π, assembler, section::TotalLagrangeSection, U, state, props, X)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)

    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    props_el = SVector{NP, eltype(props)}(@views props[:, e])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])

    R_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])
      Π_q, R_q = energy_and_internal_force(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
      Π[1] = Π[1] + Π_q
      R_el = R_el + R_q
    end

    # assemble into global matrix
    assemble!(assembler, R_el, conn)
  end
end


function energy_internal_force_and_stiffness!(Π, assembler, section::TotalLagrangeSection, U, state, props, X, block_id)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)

    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    props_el = SVector{NP, eltype(props)}(@views props[:, e])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])

    R_el = zeros(SVector{NDOF, eltype(U)})
    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])
      Π_q, R_q, K_q = energy_internal_force_and_stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
      Π[1] = Π[1] + Π_q
      R_el = R_el + R_q
      K_el = K_el + K_q
    end

    # assemble into global matrix
    assemble!(assembler, R_el, conn)
    assemble!(assembler, K_el, block_id, e)
  end
end

function energy_internal_force_and_stiffness_action!(Π, R, Hv, section::TotalLagrangeSection, U, state, props, X, V)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NS = ConstitutiveModels.num_state_vars(section.model)
  NP = ConstitutiveModels.num_properties(section.model)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)

    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    props_el = SVector{NP, eltype(props)}(@views props[:, e])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    V_el = SVector{ND * NN, eltype(V)}(@views V[conn])

    R_el = zeros(SVector{NDOF, eltype(U)})
    Hv_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])
      Π_q, R_q, H_q = energy_internal_force_and_stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
      Π[1] = Π[1] + Π_q
      R_el = R_el + R_q
      Hv_el = Hv_el + H_q * V_el
    end

    # assemble into global matrix
    assemble!(R, Hv, R_el, Hv_el, conn)
  end
end

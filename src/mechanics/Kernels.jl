# hot loop over quadrature points for a single section
function energy(section::TotalLagrangeSection, U, state, props, X, ::NoBackend)
  W = 0.0

  for e in 1:num_elements(section)
    _, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)

      W = W + energy(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
    end
  end

  return W
end

function internal_force!(f, section::TotalLagrangeSection, U, state, props, X, ::NoBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{ND * NN, eltype(f)})
    
    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)
      f_q = internal_force(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
      f_el = f_el + f_q
    end

    # assemble into global field
    assemble!(f, f_el, conn)
  end
end

function stiffness!(assembler, section::TotalLagrangeSection, U, state, props, X, block_id, ::NoBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    _, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)

      K_el = K_el + stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
    end

    # assemble into global matrix
    assemble!(assembler, K_el, block_id, e)
  end
end

function stiffness_action!(Kv, section::TotalLagrangeSection, U, state, props, X, V, ::NoBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el, V_el = unpack_element(section, U, props, X, V, e)

    Kv_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)

      K_q = stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
      Kv_el = Kv_el + K_q * V_el
    end

    # assemble into global matrix
    # assemble!(assembler, K_el, block_id, e)
    assemble!(Kv, Kv_el, conn)
  end
end

function energy_and_internal_force!(Π, f, section::TotalLagrangeSection, U, state, props, X, ::NoBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)
      Π_q, f_q = energy_and_internal_force(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
      Π[1] = Π[1] + Π_q
      f_el = f_el + f_q
    end

    # assemble into global matrix
    assemble!(f, f_el, conn)
  end
end

function internal_force_and_stiffness!(f, assembler, section::TotalLagrangeSection, U, state, props, X, block_id, ::NoBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{NDOF, eltype(f)})
    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)
      f_q, K_q = internal_force_and_stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
      f_el = f_el + f_q
      K_el = K_el + K_q
    end

    # assemble into global fields
    assemble!(f, f_el, conn)
    assemble!(assembler, K_el, block_id, e)
  end
end

function energy_internal_force_and_stiffness!(Π, f, assembler, section::TotalLagrangeSection, U, state, props, X, block_id, ::NoBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{NDOF, eltype(f)})
    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)
      Π_q, f_q, K_q = energy_internal_force_and_stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
      Π[1] = Π[1] + Π_q
      f_el = f_el + f_q
      K_el = K_el + K_q
    end

    # assemble into global fields
    assemble!(f, f_el, conn)
    assemble!(assembler, K_el, block_id, e)
  end
end

function energy_internal_force_and_stiffness_action!(Π, f, Hv, section::TotalLagrangeSection, U, state, props, X, V, ::NoBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el, V_el = unpack_element(section, U, props, X, V, e)

    f_el = zeros(SVector{NDOF, eltype(U)})
    Hv_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)
      Π_q, f_q, H_q = energy_internal_force_and_stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
      Π[1] = Π[1] + Π_q
      f_el = f_el + f_q
      Hv_el = Hv_el + H_q * V_el
    end

    # assemble into global fields
    assemble!(f, Hv, f_el, Hv_el, conn)
  end
end

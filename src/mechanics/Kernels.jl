# hot loop over quadrature points for a single section
function energy!(Πs, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, ::Backend)
  for e in 1:num_elements(section)
    _, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state_old, q, e)

      W_q, state_new_q = energy(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      Πs[q, e] = W_q
      state_new[:, q, e] .= state_new_q
    end
  end
  return nothing
end

function internal_force!(f, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, ::Backend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{ND * NN, eltype(f)})
    
    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state_old, q, e)
      f_q, state_new_q = internal_force(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      f_el = f_el + f_q
      state_new[:, q, e] .= state_new_q
    end

    # assemble into global field
    assemble!(f, f_el, conn)
  end
end

function stiffness!(assembler, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, block_id, ::Backend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    _, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state_old, q, e)

      K_q, state_new_q = stiffness(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      K_el = K_el + K_q
      state_new[:, q, e] .= state_new_q
    end

    # assemble into global matrix
    assemble!(assembler, K_el, block_id, e)
  end
end

function stiffness_action!(Kv, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, V, ::Backend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el, V_el = unpack_element(section, U, props, X, V, e)

    Kv_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state_old, q, e)

      K_q, state_new_q = stiffness(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      Kv_el = Kv_el + K_q * V_el
      state_new[:, q, e] .= state_new_q
    end

    # assemble into global matrix
    # assemble!(assembler, K_el, block_id, e)
    assemble!(Kv, Kv_el, conn)
  end
end

function energy_and_internal_force!(Πs, f, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, ::Backend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state_old, q, e)
      W_q, f_q, state_new_q = energy_and_internal_force(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      Πs[q, e] = W_q
      f_el = f_el + f_q
      state_new[:, q, e] .= state_new_q
    end

    # assemble into global matrix
    assemble!(f, f_el, conn)
  end
end

function internal_force_and_stiffness!(f, assembler, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, block_id, ::Backend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{NDOF, eltype(f)})
    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state_old, q, e)
      f_q, K_q, state_new_q = internal_force_and_stiffness(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      f_el = f_el + f_q
      K_el = K_el + K_q
      state_new[:, q, e] .= state_new_q
    end

    # assemble into global fields
    assemble!(f, f_el, conn)
    assemble!(assembler, K_el, block_id, e)
  end
end

function energy_internal_force_and_stiffness!(Πs, f, assembler, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, block_id, ::Backend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    f_el = zeros(SVector{NDOF, eltype(f)})
    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state_old, q, e)
      W_q, f_q, K_q, state_new_q = energy_internal_force_and_stiffness(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      Πs[q, e] = W_q
      f_el = f_el + f_q
      K_el = K_el + K_q
      state_new[:, q, e] .= state_new_q
    end

    # assemble into global fields
    assemble!(f, f_el, conn)
    assemble!(assembler, K_el, block_id, e)
  end
end

function energy_internal_force_and_stiffness_action!(Πs, f, Hv, state_new, section::TotalLagrangeSection, Δt, X, U, props, state_old, V, ::Backend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NDOF = ND * NN

  for e in 1:num_elements(section)
    conn, U_el, props_el, X_el, V_el = unpack_element(section, U, props, X, V, e)

    f_el = zeros(SVector{NDOF, eltype(U)})
    Hv_el = zeros(SVector{NDOF, eltype(U)})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_old_q = unpack_state(section, state, q, e)
      W_q, f_q, H_q, state_new_q = energy_internal_force_and_stiffness(
        section.model, section.formulation,
        Δt, X_el, U_el, props_el, state_old_q,
        interpolants(section, q)...
      )
      Πs[q, e] = W_q
      f_el = f_el + f_q
      Hv_el = Hv_el + H_q * V_el
      state_new[:, q, e] .= state_new_q
    end

    # assemble into global fields
    assemble!(f, Hv, f_el, Hv_el, conn)
  end
end

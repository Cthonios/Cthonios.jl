function energy!(Πs, section::TotalLagrangeSection, U, state, props, X, ::ThreadsBackend)
  Threads.@threads for e in 1:num_elements(section)
    _, U_el, props_el, X_el = unpack_element(section, U, props, X, e)

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = unpack_state(section, state, q, e)

      W_q = energy(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        interpolants(section, q)...
      )
      Πs[q, e] = W_q
    end
  end
  return nothing
end

function internal_force!(f, section::TotalLagrangeSection, U, state, props, X, ::ThreadsBackend)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)

  lck = ReentrantLock()
  Threads.@threads for e in 1:num_elements(section)
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
    Threads.lock(lck) do
      assemble!(f, f_el, conn)
    end
  end
end
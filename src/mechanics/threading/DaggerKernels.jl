function energy!(Πs, section::TotalLagrangeSection, U, state, props, X, ::DaggerBackend)

  # get sizes
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NP = ConstitutiveModels.num_properties(section.model)

  @sync for e in 1:num_elements(section)
    # out = Dagger.@spawn unpack_element(section, U, props, X, e)
    # #  _, U_el, props_el, X_el = out
    # @show out
    # @assert false

    # setup arrays
    conn = dof_connectivity(section, e)
    U_el = @views SMatrix{ND, NN, eltype(U), ND * NN}(U[conn])
    props_el = @views SVector{NP, eltype(props)}(props[:, e])
    X_el = @views SMatrix{ND, NN, eltype(X), ND * NN}(X[conn])

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
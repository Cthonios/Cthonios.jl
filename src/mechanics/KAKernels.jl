@kernel function energy_kernel!(
  Πs, 
  section::TotalLagrangeSection,
  U, state, props, X,
)
  # indices
  q, e = @index(Global, NTuple)

  # unpack arrays
  _, U_el, state_q, props_el, X_el = unpack(section, U, state, props, X, q, e)

  # get interpolants
  N, ∇N_ξ, w = interpolants(section, q)

  # run routine
  Π_q = energy(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    N, ∇N_ξ, w
  )

  # update arrays
  Πs[q, e] = Π_q

  # return nothing # issue in kernel abstractions
end

@kernel function internal_force_kernel!(
  f, 
  section::TotalLagrangeSection,
  U, state, props, X, block_id
)
  # indices
  q, e = @index(Global, NTuple)

  # unpack arrays
  conn, U_el, state_q, props_el, X_el = unpack(section, U, state, props, X, q, e)

  # get interpolants
  N, ∇N_ξ, w = interpolants(section, q)

  # run routine
  f_q = internal_force(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    N, ∇N_ξ, w
  )

  # update arrays
  FiniteElementContainers.assemble_atomic!(f, f_q, conn)
end

@kernel function stiffness_kernel!(
  assembler::StaticAssembler, 
  section::TotalLagrangeSection,
  U, state, props, X, block_id
)
  # indices
  q, e = @index(Global, NTuple)

  # unpack arrays
  _, U_el, state_q, props_el, X_el = unpack(section, U, state, props, X, q, e)

  # get interpolants
  N, ∇N_ξ, w = interpolants(section, q)

  # run routine
  K_q = stiffness(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    N, ∇N_ξ, w
  )

  # update arrays
  FiniteElementContainers.assemble_atomic!(assembler, K_q, block_id, e)
end

@kernel function internal_force_and_stiffness_kernel!(
  f,
  assembler::StaticAssembler, 
  section::TotalLagrangeSection,
  U, state, props, X, block_id
)
  # indices
  q, e = @index(Global, NTuple)

  # unpack arrays
  conn, U_el, state_q, props_el, X_el = unpack(section, U, state, props, X, q, e)

  # get interpolants
  N, ∇N_ξ, w = interpolants(section, q)

  # run routine
  f_q, K_q = internal_force_and_stiffness(
    section.model, section.formulation,
    U_el, state_q, props_el, X_el,
    N, ∇N_ξ, w
  )

  # update arrays
  FiniteElementContainers.assemble_atomic!(f, f_q, conn)
  FiniteElementContainers.assemble_atomic!(assembler, K_q, block_id, e)
end


# @kernel function internal_force_and_stiffness_kernel!(assembler, section::TotalLagrangeSection, U, state, props, X, block_id)
#   q, e = @index(Global, NTuple)

#   # stuff for static arrays
#   ND = FiniteElementContainers.num_dimensions(section)
#   NN = FiniteElementContainers.num_nodes_per_element(section)
#   NS = ConstitutiveModels.num_state_vars(section.model)
#   NP = ConstitutiveModels.num_properties(section.model)

#   # unpack arrays
#   conn = dof_connectivity(section, e)

#   U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
#   props_el = SVector{NP, eltype(props)}(@views props[:, e])
#   X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
#   state_q = SVector{NS, eltype(props)}(@views state[:, q, e])

#   # run routine
#   R_q, K_q = internal_force_and_stiffness(
#     section.model, section.formulation,
#     U_el, state_q, props_el, X_el,
#     shape_function_values(section, q),
#     shape_function_gradients(section, q),
#     quadrature_weights(section, q)
#   )

#   # update arrays
#   FiniteElementContainers.assemble_atomic!(assembler, R_q, conn)
#   FiniteElementContainers.assemble_atomic!(assembler, K_q, block_id, e)

#   return nothing
# end
# energy methods below

# lazy method for developers in the REPL
function energy(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return energy(domain, Uu, state, props, X)
end    

function internal_force(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return internal_force(domain, Uu, state, props, X)
end   

function stiffness(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return stiffness(domain, Uu, state, props, X)
end

function internal_force_and_stiffness(domain::QuasiStaticDomain)
  # unpack stuff
  Uu = domain.domain_cache.Uu
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  X = domain.coords

  return internal_force_and_stiffness(domain, Uu, state, props, X)
end

# out of place methods with differentiable stuff exposed
function energy(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  Π = zeros(eltype(Uu), 1)
  energy!(Π, U, domain, Uu, state, props, X)
  return Π
end

function internal_force(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  R = create_fields(domain)
  internal_force!(R, U, domain, Uu, state, props, X)
  return R
end

function stiffness(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  assembler = domain.assembler
  stiffness!(assembler, U, domain, Uu, state, props, X)
  K = SparseArrays.sparse!(assembler)
  return 0.5 * (K + K')
end

function internal_force_and_stiffness(domain::QuasiStaticDomain, Uu::V, state, props, X) where V <: AbstractVector
  U = create_fields(domain)
  assembler = domain.assembler
  internal_force_and_stiffness!(assembler, U, domain, Uu, state, props, X)
  # R = @views assembler.residuals[domain.dof.unknown_dofs]
  K = SparseArrays.sparse!(assembler)
  return assembler.residuals, 0.5 * (K + K')
end

# top level in place method that loops over sections
function energy!(Π, U, domain::QuasiStaticDomain, Uu, state, props, X)
  update_fields!(U, domain, X, Uu)
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    Π[1] = Π[1] + energy(section, U, state_temp, props_temp, X)
  end 
  return nothing
end 

function internal_force!(R, U, domain::QuasiStaticDomain, Uu, state, props, X)
  update_fields!(U, domain, X, Uu)
  R .= zero(eltype(R))
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    internal_force!(R, section, U, state_temp, props_temp, X)
  end 
  return nothing
end 

function stiffness!(assembler, U, domain::QuasiStaticDomain, Uu, state, props, X)
  update_fields!(U, domain, X, Uu)
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    stiffness!(assembler, section, U, state_temp, props_temp, X, block_count)
    block_count = block_count + 1
  end 
  return nothing
end 



# function internal_force_and_stiffness_old!(assembler, U, domain::QuasiStaticDomain, Uu, state, props, X)
#   update_fields!(U, domain, X, Uu)
#   assembler.residuals .= zero(eltype(assembler.residuals))
#   block_count = 1 # needed for the in place assembler from FiniteElementContainers
#   for (name, section) in pairs(domain.sections)
#     state_temp = @views state[name]
#     props_temp = @views props[name]
#     internal_force_and_stiffness!(assembler, section, U, state_temp, props_temp, X, block_count)
#     block_count = block_count + 1
#   end 
#   return nothing
# end 

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

function internal_force_and_stiffness!(assembler, section::TotalLagrangeSection, U, state, props, X, block_id)
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
    R_el = zeros(SVector{ND * NN, eltype(assembler.residuals)})
    K_el = zeros(SMatrix{NDOF, NDOF, eltype(U), NDOF * NDOF})

    for q in 1:FiniteElementContainers.num_q_points(section)
      state_q = SVector{NS, eltype(props)}(@views state[:, q, e])
      R_q, K_q = internal_force_and_stiffness(
        section.model, section.formulation,
        U_el, state_q, props_el, X_el,
        shape_function_values(section, q),
        shape_function_gradients(section, q),
        quadrature_weights(section, q)
      )
      R_el = R_el + R_q
      K_el = K_el + K_q
    end

    # assemble into global matrix
    assemble!(assembler, R_el, conn)
    assemble!(assembler, K_el, block_id, e)
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

function internal_force_and_stiffness!(assembler, U, domain::QuasiStaticDomain, Uu, state, props, X)
  backend = CPU()
  kernel! = internal_force_and_stiffness_kernel!(backend)

  update_fields!(U, domain, X, Uu)
  assembler.residuals .= zero(eltype(assembler.residuals))
  block_count = 1 # needed for the in place assembler from FiniteElementContainers
  for (name, section) in pairs(domain.sections)
    state_temp = @views state[name]
    props_temp = @views props[name]
    # internal_force_and_stiffness!(assembler, section, U, state_temp, props_temp, X, block_count)
    kernel!(assembler, section, U, state_temp, props_temp, X, block_count, ndrange=(FiniteElementContainers.num_q_points(section), num_elements(section)))
    block_count = block_count + 1
  end 
  return nothing
end 

# quadrature level methods, the actual mechanics is here
function energy(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  U_el, state_q, props_el, X_el, 
  N, ∇N_ξ, w
)

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  ψ_q = ConstitutiveModels.energy(model, props_el, F_q, state_q)

  return w * det(J) * ψ_q
end

function internal_force(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  U_el, state_q, props_el, X_el, 
  N, ∇N_ξ, w
)

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  # ψ_q = ConstitutiveModels.energy(model, props_el, F_q, state_q)
  P_q = ConstitutiveModels.pk1_stress(model, props_el, F_q, state_q)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return w * det(J) * G * P_v
end

function stiffness(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  U_el, state_q, props_el, X_el, 
  N, ∇N_ξ, w
)

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  A_q = Tensors.gradient(x -> ConstitutiveModels.pk1_stress(model, props_el, x, state_q), F_q)
  A_v = FiniteElementContainers.extract_stiffness(formulation, A_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return w * det(J) * G * A_v * G'
end

function internal_force_and_stiffness(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  U_el, state_q, props_el, X_el, 
  N, ∇N_ξ, w
)

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  A_q, P_q = Tensors.gradient(x -> ConstitutiveModels.pk1_stress(model, props_el, x, state_q), F_q, :all)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q)
  A_v = FiniteElementContainers.extract_stiffness(formulation, A_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)

  return w * det(J) * G * P_v, w * det(J) * G * A_v * G'
end

function gradient end
function gradient! end
function residual end
function residual! end
function residual_dot_v! end
function hvp end
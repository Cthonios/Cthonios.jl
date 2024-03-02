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

function energy_and_internal_force(
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
  # TODO eventually have a method in ConstitutiveModels that does this
  # TODO and intelligently uses AD where we don't have analytic implementations
  P_q, ψ_q = Tensors.gradient(x -> ConstitutiveModels.energy(model, props_el, x, state_q), F_q, :all)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q)
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)

  return w * det(J) * ψ_q, w * det(J) * G * P_v
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

function energy_internal_force_and_stiffness(
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
  # TODO eventually have a method in ConstitutiveModels that does this
  # TODO and intelligently uses AD where we don't have analytic implementations
  A_q, P_q, ψ_q = Tensors.hessian(x -> ConstitutiveModels.energy(model, props_el, x, state_q), F_q, :all)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q)
  A_v = FiniteElementContainers.extract_stiffness(formulation, A_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)

  return w * det(J) * ψ_q, w * det(J) * G * P_v, w * det(J) * G * A_v * G'
end

function energy_gradient end

include("Utils.jl")
include("Kernels.jl")
# include("DaggerKernels.jl")
# include("KAKernels.jl")
# include("ThreadsKernels.jl")
include("SectionIterators.jl")
include("DomainWrappers.jl")
include("Sensitivities.jl")

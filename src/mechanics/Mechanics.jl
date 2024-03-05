# quadrature level methods, the actual mechanics is here



function energy(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  Δt, X_el, U_el, props_el, state_old_q,
  N, ∇N_ξ, w
)

  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  ψ_q, state_new_q = ConstitutiveModels.helmholtz_free_energy(model, props_el, Δt, F_q, θ_q, state_old_q)

  return w * det(J) * ψ_q, state_new_q
end

function internal_force(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  Δt, X_el, U_el, props_el, state_old_q,
  N, ∇N_ξ, w
)

  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  # ψ_q = ConstitutiveModels.energy(model, props_el, F_q, state_q)
  P_q, state_new_q = ConstitutiveModels.pk1_stress(model, props_el, Δt, F_q, θ_q, state_old_q)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return w * det(J) * G * P_v, state_new_q
end

function stiffness(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  Δt, X_el, U_el, props_el, state_old_q,
  N, ∇N_ξ, w
)

  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  # A_q, state_new_q = Tensors.gradient(x -> ConstitutiveModels.pk1_stress(model, props_el, Δt, F_q, θ_q, state_old_q), F_q)
  A_q, state_new_q = ConstitutiveModels.material_tangent(model, props_el, Δt, F_q, θ_q, state_old_q)
  A_v = FiniteElementContainers.extract_stiffness(formulation, A_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return w * det(J) * G * A_v * G', state_new_q
end

function energy_and_internal_force(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  Δt, X_el, U_el, props_el, state_old_q,
  N, ∇N_ξ, w
)

  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0

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
  # P_q, ψ_q = Tensors.gradient(x -> ConstitutiveModels.energy(model, props_el, Δt, F_q, θ_q, state_old_q), F_q, :all)
  ψ_q, P_q, state_new_q = ConstitutiveModels.helmholtz_free_energy_and_pk1_stress(model, props_el, Δt, F_q, θ_q, state_old_q)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q)
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)

  return w * det(J) * ψ_q, w * det(J) * G * P_v, state_new_q
end

function internal_force_and_stiffness(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  Δt, X_el, U_el, props_el, state_old_q,
  N, ∇N_ξ, w
)

  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0

  # map shape functions
  J    = (X_el * ∇N_ξ)'
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)

  # constitutive
  # A_q, P_q = Tensors.gradient(x -> ConstitutiveModels.pk1_stress(model, props_el, x, state_q), F_q, :all)
  P_q, A_q, state_new_q = ConstitutiveModels.pk1_stress_and_material_tangent(model, props_el, Δt, F_q, θ_q, state_old_q)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q)
  A_v = FiniteElementContainers.extract_stiffness(formulation, A_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)

  return w * det(J) * G * P_v, w * det(J) * G * A_v * G', state_new_q
end

function energy_internal_force_and_stiffness(
  model::ConstitutiveModels.MechanicalModel, 
  formulation::AbstractMechanicsFormulation,
  Δt, X_el, U_el, props_el, state_old_q,
  N, ∇N_ξ, w
)

  # TODO temp hardcoded time step and temperature
  # Hook up stream with everything elase
  θ_q = 0.0

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
  # A_q, P_q, ψ_q = Tensors.hessian(x -> ConstitutiveModels.energy(model, props_el, x, state_q), F_q, :all)
  ψ_q, P_q, A_q, state_new_q = ConstitutiveModels.mechanical_state(model, props_el, Δt, F_q, θ_q, state_old_q)
  P_v = FiniteElementContainers.extract_stress(formulation, P_q)
  A_v = FiniteElementContainers.extract_stiffness(formulation, A_q) 
  G = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)

  return w * det(J) * ψ_q, w * det(J) * G * P_v, w * det(J) * G * A_v * G', state_new_q
end

function energy_gradient end
function energy_hvp end

include("Utils.jl")
include("SectionIterators.jl")
include("DomainWrappers.jl")
# include("Sensitivities.jl")

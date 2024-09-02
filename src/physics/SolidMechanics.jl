"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct SolidMechanics{
  NF, NP, NS,
  Mat <: ConstitutiveModel{NP, NS}, 
  Form <: AbstractMechanicsFormulation{NF}
} <: AbstractPhysics{NF, NP, NS}
  material_model::Mat
  formulation::Form
end

# function discrete_gradient(physics::SolidMechanics, ∇N_X)
#   return FiniteElementContainers.discrete_gradient(
#     physics.formulation, ∇N_X
#   )
# end

# function extract_stiffness(physics::SolidMechanics, A)
#   return FiniteElementContainers.extract_stiffness(
#     physics.formulation, A
#   )
# end 

# function extract_stress(physics::SolidMechanics, P)
#   return FiniteElementContainers.extract_stress(
#     physics.formulation, P
#   )
# end

# function modify_field_gradients(physics::SolidMechanics, ∇u)
#   return FiniteElementContainers.modify_field_gradients(
#     physics.formulation, ∇u
#   )
# end

function energy(physics::SolidMechanics, cell, u_el)
  @unpack X_q, N, ∇N_X, JxW = cell
  ∇u_q = u_el * ∇N_X
  ∇u_q = modify_field_gradients(physics/formulation, ∇u_q)
  F_q = ∇u_q + one(∇u_q)

  # hardcoded for now
  props = SVector{2, Float64}((0.833, 0.3846))
  dt = 0.0
  θ = 0.0
  Q = SVector{0, Float64}()

  ψ, Q = ConstitutiveModels.helmholtz_free_energy(
    physics.material_model, props, dt, F_q, θ, Q
  )
  return JxW * Ψ
end

function gradient(physics::SolidMechanics, cell, u_el)
  @unpack X_q, N, ∇N_X, JxW = cell
  ∇u_q = u_el * ∇N_X
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  F_q = ∇u_q + one(∇u_q)

  # hardcoded for now
  props = SVector{2, Float64}((0.833, 0.3846))
  dt = 0.0
  θ = 0.0
  Q = SVector{0, Float64}()

  # constitutive
  P_q, state_new_q = ConstitutiveModels.pk1_stress(
    physics.material_model, props, dt, F_q, θ, Q
  )
  P_v = extract_stress(physics.formulation, P_q) 
  G = discrete_gradient(physics.formulation, ∇N_X)
  return JxW * G * P_v#, state_new_q
end

function hessian(physics::SolidMechanics, cell, u_el)
  @unpack X_q, N, ∇N_X, JxW = cell
  ∇u_q = u_el * ∇N_X
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  F_q = ∇u_q + one(∇u_q)

  # hardcoded for now
  props = SVector{2, Float64}((0.833, 0.3846))
  dt = 0.0
  θ = 0.0
  Q = SVector{0, Float64}()

  # constitutive
  A_q, state_new_q = ConstitutiveModels.material_tangent(
    physics.material_model, props, dt, F_q, θ, Q
  )
  A_v = extract_stiffness(physics.formulation, A_q) 
  G = discrete_gradient(physics.formulation, ∇N_X)
  return JxW * G * A_v * G'#, state_new_q
end
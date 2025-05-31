"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct StressDivergence{
  NF, NP, NS,
  Mat <: ConstitutiveModel{NP, NS}, 
  Form <: AbstractMechanicsFormulation{NF}
} <: AbstractPhysics{NF, NP, NS}
  material_model::Mat
  formulation::Form
end

init_properties(physics::StressDivergence, props) = ConstitutiveModels.initialize_props(physics.material_model, props)
init_state(physics::StressDivergence) = ConstitutiveModels.initialize_state(physics.material_model)

"""
$(TYPEDSIGNATURES)
Energy method at the quadrature level for
Lagrangian solid mechanics. This equivalent to
the quadrature point calculation needed for the 
following integral
``
\\Pi = \\int_\\Omega\\psi\\left(\\mathbf{F}\\right)d\\Omega
``
"""
function energy(physics::StressDivergence, u::T, ∇u, X, t, Z, props) where T <: AbstractArray
  F = ∇u + one(∇u)
  dt = time_step(t)

  # hardcoded for now
  θ = 0.0

  ψ, Q = ConstitutiveModels.helmholtz_free_energy(
    physics.material_model, props, dt, F, θ, Z
  )
  return ψ
end

"""
$(TYPEDSIGNATURES)
Gradient method at the quadrature level for
Lagrangian solid mechanics. This equivalent to
the quadrature point calculation needed for the 
following integral
``
\\mathbf{f} = \\int_\\Omega\\mathbf{P}:\\delta\\mathbf{F}d\\Omega
``
"""
function gradient(physics::StressDivergence, u, ∇u, v, ∇v, X, t, Z, props)
  F = ∇u + one(∇u)
  dt = time_step(t)

  # hardcoded for now
  θ = 0.0

  # constitutive
  P, Q = ConstitutiveModels.pk1_stress(
    physics.material_model, props, dt, F, θ, Z
  )
  P = extract_stress(physics.formulation, P)
  return ∇v * P
end

"""
$(TYPEDSIGNATURES)
"""
function hessian(physics::StressDivergence, u, ∇u, v, ∇v, X, t, Z, props)
  F = ∇u + one(∇u)
  dt = time_step(t)

  # hardcoded for now
  θ = 0.0

  # constitutive
  A, Q = ConstitutiveModels.material_tangent(
    physics.material_model, props, dt, F, θ, Z
  )
  A = extract_stiffness(physics.formulation, A)
  return ∇v * A * ∇v'
end

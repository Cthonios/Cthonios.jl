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

function SolidMechanics(inputs::Dict{Symbol, Any})
  material_inputs = inputs[:material]
  model_name = material_inputs[:type]
  model = eval(Meta.parse(model_name))(material_inputs[:properties])
  formulation_inputs = inputs[:formulation]
  formulation = eval(Symbol(formulation_inputs[:type]))()
  return SolidMechanics(model, formulation)
end

init_properties(physics::SolidMechanics, props) = ConstitutiveModels.initialize_props(physics.material_model, props)

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
function energy(physics::SolidMechanics, u::T, ∇u, X, props) where T <: AbstractArray
  F = ∇u + one(∇u)

  # hardcoded for now
  dt = 0.0
  θ = 0.0
  Q = SVector{0, Float64}()

  ψ, Q = ConstitutiveModels.helmholtz_free_energy(
    physics.material_model, props, dt, F, θ, Q
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
function gradient(physics::SolidMechanics, u, ∇u, v, ∇v, X, props)
  F = ∇u + one(∇u)

  # hardcoded for now
  dt = 0.0
  θ = 0.0
  Q = SVector{0, Float64}()

  # constitutive
  P, Q = ConstitutiveModels.pk1_stress(
    physics.material_model, props, dt, F, θ, Q
  )
  P = extract_stress(physics.formulation, P)
  return ∇v * P
end

"""
$(TYPEDSIGNATURES)
"""
function hessian(physics::SolidMechanics, u, ∇u, v, ∇v, X, props)
  F = ∇u + one(∇u)

  # hardcoded for now
  dt = 0.0
  θ = 0.0
  Q = SVector{0, Float64}()

  # constitutive
  A, Q = ConstitutiveModels.material_tangent(
    physics.material_model, props, dt, F, θ, Q
  )
  A = extract_stiffness(physics.formulation, A)
  return ∇v * A * ∇v'
end

function neumann_bc_energy(physics::SolidMechanics)

end

function neumann_bc_gradient()

end
"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct SolidMechanics{
  NF, NP, NS, Mat, Form,
  P <: StressDivergence{NF, NP, NS, Mat, Form}
} <: AbstractPhysics{NF, NP, NS}
  # TODO redundant fields below since
  # material model and formulation are both in
  # stress divergence
  material_model::Mat
  formulation::Form
  stress_divergence::P
end

# script constructor
function SolidMechanics(model::ConstitutiveModel, form::AbstractMechanicsFormulation)
  stress_divergence = StressDivergence(model, form)
  return SolidMechanics(model, form, stress_divergence)
end

# input file constructor
function SolidMechanics(inputs::Dict{Symbol, Any})
  material_inputs = inputs[:material]
  model_name = material_inputs[:type]
  model = eval(Meta.parse(model_name))(material_inputs[:properties])
  formulation_inputs = inputs[:formulation]
  formulation = eval(Symbol(formulation_inputs[:type]))()
  stress_divergence = StressDivergence(model, formulation)
  return SolidMechanics(formulation, model, stress_divergence)
end

init_properties(physics::SolidMechanics, props) = init_properties(physics.stress_divergence, props)
init_state(physics::SolidMechanics) = init_state(physics.stress_divergence)

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
function energy(physics::SolidMechanics, u::T, ∇u, X, t, Z, props) where T <: AbstractArray
  internal_energy = energy(physics.stress_divergence, u, ∇u, X, t, Z, props)
  return internal_energy
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
function gradient(physics::SolidMechanics, u, ∇u, v, ∇v, X, t, Z, props)
  internal_force = gradient(physics.stress_divergence, u, ∇u, v, ∇v, X, t, Z, props)
  return internal_force
end

"""
$(TYPEDSIGNATURES)
"""
function hessian(physics::SolidMechanics, u, ∇u, v, ∇v, X, t, Z, props)
  stiffness = hessian(physics.stress_divergence, u, ∇u, v, ∇v, X, t, Z, props)
  return stiffness
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct SolidDynamics{
  NF, NP, NS, Mat, Form,
  P1 <: StressDivergence{NF, NP, NS, Mat, Form},
  P2 <: Dynamics{NF}
} <: AbstractPhysics{NF, NP, NS}
  # TODO redundant fields below since
  # material model and formulation are both in
  # stress divergence
  material_model::Mat
  formulation::Form
  stress_divergence::P1
  dynamics::P2
end

# script constructor
function SolidDynamics(model::ConstitutiveModel, form::AbstractMechanicsFormulation)
  stress_divergence = StressDivergence(model, form)
  dynamics = Dynamics{FiniteElementContainers.num_dimensions(form)}()
  return SolidDynamics(model, form, stress_divergence, dynamics)
end

# input file constructor
function SolidDynamics(inputs::Dict{Symbol, Any})
  material_inputs = inputs[:material]
  model_name = material_inputs[:type]
  model = eval(Meta.parse(model_name))(material_inputs[:properties])
  formulation_inputs = inputs[:formulation]
  formulation = eval(Symbol(formulation_inputs[:type]))()
  stress_divergence = StressDivergence(model, formulation)
  dynamics = Dynamics{FiniteElementContainers.num_dimensions(formulation)}()
  return SolidDynamics(formulation, model, stress_divergence, dynamics)
end

init_properties(physics::SolidDynamics, props) = init_properties(physics.stress_divergence, props)
init_state(physics::SolidDynamics) = init_state(physics.stress_divergence)

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
function energy(physics::SolidDynamics, u::T, ∇u, X, t, Z, props) where T <: AbstractArray
  internal_energy = energy(physics.stress_divergence, u, ∇u, X, t, Z, props)
  kinetic_energy = energy(physics.dynamics, u, ∇u, X, t, Z, props)
  return internal_energy + kinetic_energy
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
function gradient(physics::SolidDynamics, u, ∇u, v, ∇v, X, t, Z, props)
  internal_force = gradient(physics.stress_divergence, u, ∇u, v, ∇v, X, t, Z, props)
  # kinetic_force = gradient(physics.dynamics, u, ∇u, v, ∇v, X, t, Z, props)
  return internal_force
end

"""
$(TYPEDSIGNATURES)
"""
function hessian(physics::SolidDynamics, u, ∇u, v, ∇v, X, t, Z, props)
  stiffness = hessian(physics.stress_divergence, u, ∇u, v, ∇v, X, t, Z, props)
  # kinetic_force = gradient(physics.dynamics, u, ∇u, v, ∇v, X, t, Z, props)
  return stiffness
end

# function kinetic_energy(physics::SolidDynamics, u, ∇u, X, t, Z, props)
#   return energy(physics.dynamics, u, ∇u, X, t, Z, props)
# end

"""
$(TYPEDSIGNATURES)
"""
function mass_matrix(physics::SolidDynamics, u, ∇u, v, ∇v, X, t, Z, props)
  masses = hessian(physics.dynamics, u, ∇u, v, ∇v, X, t, Z, props)
  # kinetic_force = gradient(physics.dynamics, u, ∇u, v, ∇v, X, t, Z, props)
  return masses
end
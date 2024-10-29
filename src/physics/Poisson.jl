# move below to sub folders maybe?
"""
$(TYPEDEF)
$(TYPEDFIELDS)
# TODO create source term kernel
"""
struct Poisson{Form, F} <: AbstractPhysics{1, 0, 0}
  formulation::Form
  laplacian::Laplacian{1}
  func::F
end

function Poisson(f)
  formulation = ScalarFormulation()
  laplacian = Laplacian{1}()
  return Poisson{typeof(formulation), typeof(f)}(formulation, laplacian, f)
end

init_properties(physics::Poisson, props) = zeros(num_properties(physics))
init_state(physics::Poisson) = zeros(num_states(physics))

"""
$(TYPEDSIGNATURES)
Energy method for Poisson equation at a quadrature point
``
\\Pi\\left[u\\right] = \\int_\\Omega \\left[\\frac{1}{2}\\|\\nabla u\\|^2 - fu\\right]d\\Omega
``
"""
function energy(physics::Poisson, u::T, ∇u, X, t, Z, props) where T <: AbstractArray
  Π = 0.5 * dot(∇u, ∇u) - physics.func(X, 0.0) * u
  return Π
end

"""
$(TYPEDSIGNATURES)
Gradient method for Poisson equation at a quadrature point
``
g\\left(u, v\\right) = \\int_\\Omega \\left[\\nabla u\\cdot\\nabla v - fv\\right]d\\Omega
``
"""
function gradient(physics::Poisson, u, ∇u, v, ∇v, X, t, Z, props)
  R = ∇u * ∇v' - v' * physics.func(X, 0.0)
  return R[:]
end

"""
$(TYPEDSIGNATURES)
Hessian method for Poisson equation at a quadrature point
``
H\\left(u, v\\right) = \\int_\\Omega \\left[\\nabla v\\cdot\\nabla v\\right]d\\Omega
``
"""
function hessian(physics::Poisson, u, ∇u, v, ∇v, X, t, Z, props)
  K = ∇v * ∇v'
  return K
end

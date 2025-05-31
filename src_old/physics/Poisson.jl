# move below to sub folders maybe?
"""
$(TYPEDEF)
$(TYPEDFIELDS)
# TODO create source term kernel
"""
struct Poisson{Form, S} <: AbstractPhysics{1, 0, 0}
  formulation::Form
  laplacian::Laplacian{1}
  source::S
end

function Poisson(f)
  formulation = ScalarFormulation()
  laplacian = Laplacian{1}()
  source = Source{1}(f)
  return Poisson(formulation, laplacian, source)
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
  Π = energy(physics.laplacian, u, ∇u, X, t, Z, props) + 
      energy(physics.source, u, ∇u, X, t, Z, props)
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
  R = gradient(physics.laplacian, u, ∇u, v, ∇v, X, t, Z, props) + 
      gradient(physics.source, u, ∇u, v, ∇v, X, t, Z, props)
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
  K = hessian(physics.laplacian, u, ∇u, v, ∇v, X, t, Z, props) + 
      hessian(physics.source, u, ∇u, v, ∇v, X, t, Z, props)
  return K
end

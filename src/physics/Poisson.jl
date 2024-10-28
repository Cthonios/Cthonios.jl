# move below to sub folders maybe?
"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Poisson{F} <: AbstractPhysics{1, 0, 0}
  func::F
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
function energy(physics::Poisson, u::T, ∇u, X, props) where T <: AbstractArray
  # u_el = reshape_field(physics, cell, u_el)
  # @unpack X_q, N, ∇N_X, JxW = cell
  # u_q = u_el * N
  # ∇u_q = u_el * ∇N_X
  Π = 0.5 * dot(∇u, ∇u) - physics.func(X, 0.0) * u
  # return cell.JxW * Π
  return Π
end

"""
$(TYPEDSIGNATURES)
Gradient method for Poisson equation at a quadrature point
``
g\\left(u, v\\right) = \\int_\\Omega \\left[\\nabla u\\cdot\\nabla v - fv\\right]d\\Omega
``
"""
function gradient(physics::Poisson, u, ∇u, v, ∇v, X, props)
  # u_el = reshape_field(physics, cell, u_el)
  # @unpack X_q, N, ∇N_X, JxW = cell
  # ∇u_q = u_el * ∇N_X
  R = ∇u * ∇v' - v' * physics.func(X, 0.0)
  # return JxW * R[:]
  return R
end

"""
$(TYPEDSIGNATURES)
Hessian method for Poisson equation at a quadrature point
``
H\\left(u, v\\right) = \\int_\\Omega \\left[\\nabla v\\cdot\\nabla v\\right]d\\Omega
``
"""
function hessian(physics::Poisson, u, ∇u, v, ∇v, X, props)
  # u_el = reshape_field(physics, cell, u_el)
  # @unpack X_q, N, ∇N_X, JxW = cell
  K = ∇v * ∇v'
  # return JxW * K_q
  return K
end

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
function energy(physics::Poisson, cell, u_el, _)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q = u_el * N
  ∇u_q = u_el * ∇N_X
  Π_q = 0.5 * dot(∇u_q, ∇u_q) - physics.func(X_q, 0.0) * u_q
  return JxW * Π_q
end

"""
$(TYPEDSIGNATURES)
Gradient method for Poisson equation at a quadrature point
``
g\\left(u, v\\right) = \\int_\\Omega \\left[\\nabla u\\cdot\\nabla v - fv\\right]d\\Omega
``
"""
function gradient(physics::Poisson, cell, u_el, _)
  @unpack X_q, N, ∇N_X, JxW = cell
  ∇u_q = u_el * ∇N_X
  R_q = ∇u_q * ∇N_X' - N' * physics.func(X_q, 0.0)
  return JxW * R_q[:]
end

"""
$(TYPEDSIGNATURES)
Hessian method for Poisson equation at a quadrature point
``
H\\left(u, v\\right) = \\int_\\Omega \\left[\\nabla v\\cdot\\nabla v\\right]d\\Omega
``
"""
function hessian(::Poisson, cell, u_el, _)
  @unpack X_q, N, ∇N_X, JxW = cell
  K_q = ∇N_X * ∇N_X'
  return JxW * K_q
end

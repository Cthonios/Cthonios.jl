"""
$(TYPEDEF)
$(TYPEDFIELDS)
Objective type for evaluating objective functions.
The functions correspond to quadrature level evaluations
of the objective function, it's gradient, and it's hessian
respectively.
```domain``` - A domain object 
```value``` - A function for the quadrature level objective function evaluation.
```gradient``` - A function for the quadrature level objective function gradient evaluation.
```hessian``` - A function for the quadrature level objective hessian evaluation.
```timer``` - A timer that's already setup.
"""
struct UnconstrainedObjective{
  D, 
  # F1, F2, F3, F4,
  V,
  # F5, F6, F7, F8, 
  S,
  T
} <: AbstractObjective
  domain::D
  # value::F1
  # gradient::F2
  # gradient_x::F3
  # hessian::F4
  volume_integral::V
  # neumann_value::F5
  # neumann_gradient::F6
  # neumann_gradient_x::F7
  # neumann_hessian::F8
  surface_integral::S
  timer::T
end

# function UnconstrainedObjective(domain::Domain, value_func::F1, neumann_value_func::F2, timer::TimerOutput) where {F1 <: Function, F2 <: Function}
#   grad_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> value_func(phys, cell, z, x, state, props, t), u)
#   grad_x_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> value_func(phys, cell, u, z, state, props, t), x)
#   hess_func(phys, cell, u, x, state, props, t) = ForwardDiff.hessian(z -> value_func(phys, cell, z, x, state, props, t), u)
#   neumann_grad_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_value_func(phys, cell, z, x, t, bc, val, r, q, f), u)
#   neumann_grad_x_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_value_func(phys, cell, u, z, t, bc, val, r, q, f), x)
#   neumann_hess_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.hessian(z -> neumann_value_func(phys, cell, z, x, t, bc, val, r, q, f), u)
#   return UnconstrainedObjective(
#     domain, 
#     value_func, grad_func, grad_x_func, hess_func, 
#     neumann_value_func, neumann_grad_func, neumann_grad_x_func, neumann_hess_func,
#     timer
#   )
# end

function UnconstrainedObjective(domain::Domain, ::typeof(energy))
  # grad_x_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> energy(phys, cell, u, z, state, props, t), x)
  neumann_grad_x_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_energy(phys, cell, u, z, t, bc, val, r, q, f), x)
  return UnconstrainedObjective(
    domain,
    # energy, gradient, grad_x_func, hessian,
    VolumeIntegral(energy),
    # neumann_energy, neumann_gradient, neumann_grad_x_func, neumann_hessian,
    SurfaceIntegral(neumann_energy),
    TimerOutput()
  )
end

function UnconstrainedObjective(inputs::Dict{Symbol, Any}, domain, timer)
  value_func = eval(Meta.parse(inputs[:value]))

  if typeof(value_func) == typeof(energy)
    return UnconstrainedObjective(domain, value_func)
  else
    throw(ErrorException("Unsupported objective function"))
  end
end

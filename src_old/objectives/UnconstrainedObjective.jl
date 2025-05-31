"""
$(TYPEDEF)
$(TYPEDFIELDS)
Objective type for evaluating objective functions.
The functions correspond to quadrature level evaluations
of the objective function, it's gradient, and it's hessian
respectively.
```domain``` - A domain object 
```integral``` - An AbstractIntegral
```timer``` - A timer that's already setup.
"""
struct UnconstrainedObjective{D <: Domain, I <: AbstractIntegral, T} <: AbstractObjective
  domain::D
  integral::I
  timer::T
end

function UnconstrainedObjective(domain::Domain, ::typeof(energy))
  # grad_x_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> energy(phys, cell, u, z, state, props, t), x)
  neumann_grad_x_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_energy(phys, cell, u, z, t, bc, val, r, q, f), x)
  return UnconstrainedObjective(
    domain,
    VariationalIntegral(
      VolumeIntegral(energy),
      SurfaceIntegral(neumann_energy),
    ),
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

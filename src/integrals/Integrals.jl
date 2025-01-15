abstract type AbstractIntegral end

function integrate! end

function gradient_x!(val, integral::AbstractIntegral, Uu, p, domain)
  integrate!(val, integral, integral.integrand_gradient_x, Uu, p, domain)
end

function gradient_u!(val, integral::AbstractIntegral, Uu, p, domain)
  integrate!(val, integral, integral.integrand_gradient_u, Uu, p, domain)
end

function hessian_u!(val, integral::AbstractIntegral, Uu, p, domain)
  integrate!(val, integral, integral.integrand_hessian_u, Uu, p, domain)
end

function hvp_u!(val, integral::AbstractIntegral, Uu, p, Vv, domain)
  integrate!(val, integral, integral.integrand_hessian_u, Uu, p, Vv, domain)
end

function value!(val, integral::AbstractIntegral, Uu, p, domain)
  integrate!(val, integral, integral.integrand, Uu, p, domain)
end

include("Utils.jl")
include("SurfaceIntegral.jl")
include("VolumeIntegral.jl")

# function Base.:+()

# end
# implementing sometihg more general

struct Integral{A, B} <: AbstractIntegral
  volume_integral::A
  surface_integral::B
end

function Integral(::typeof(energy))
  # create volume integral
  grad_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> energy(phys, cell, z, x, state, props, t), u)
  grad_x_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> energy(phys, cell, u, z, state, props, t), x)
  hess_func(phys, cell, u, x, state, props, t) = ForwardDiff.hessian(z -> energy(phys, cell, z, x, state, props, t), u)


  neumann_grad_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_energy(phys, cell, z, x, t, bc, val, r, q, f), u)
  neumann_grad_x_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_energy(phys, cell, u, z, t, bc, val, r, q, f), x)
  neumann_hess_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.hessian(z -> neumann_energy(phys, cell, z, x, t, bc, val, r, q, f), u)
end

function integrate!(val, integral::Integral, Uu, p, domain)
  update_field_unknowns!(current_solution(p), domain, Uu)
  if integral.volume_integral !== nothing
    integrate!(val, integral.volume_integral, Uu, p, domain)
  end

  if integral.surface_integral !== nothing
    integrate!(val, integral.surface_integral, Uu, p, domain)
  end
end

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

struct VariationalIntegral{A, B} <: AbstractIntegral
  volume_integral::A
  surface_integral::B
end

function integrate!(val, integral::VariationalIntegral, Uu, p, domain)
  update_field_unknowns!(current_solution(p), domain, Uu)
  if integral.volume_integral !== nothing
    integrate!(val, integral.volume_integral, Uu, p, domain)
  end

  if integral.surface_integral !== nothing
    integrate!(val, integral.surface_integral, Uu, p, domain)
  end
end

function gradient_x!(val, integral::VariationalIntegral, Uu, p, domain)
  integrate!(val, integral.volume_integral, integral.volume_integral.integrand_gradient_x, Uu, p, domain)
  integrate!(val, integral.surface_integral, integral.surface_integral.integrand_gradient_x, Uu, p, domain)
end

function gradient_u!(val, integral::VariationalIntegral, Uu, p, domain)
  integrate!(val, integral.volume_integral, integral.volume_integral.integrand_gradient_u, Uu, p, domain)
  integrate!(val, integral.surface_integral, integral.surface_integral.integrand_gradient_u, Uu, p, domain)
end

function hessian_u!(val, integral::VariationalIntegral, Uu, p, domain)
  integrate!(val, integral.volume_integral, integral.volume_integral.integrand_hessian_u, Uu, p, domain)
  # integrate!(val, integral.surface_integral, integral.surface_integral.integrand_hessian_u, Uu, p, domain)
end

function hvp_u!(val, integral::VariationalIntegral, Uu, p, Vv, domain)
  integrate!(val, integral.volume_integral, integral.volume_integral.integrand_hessian_u, Uu, p, Vv, domain)
  # integrate!(val, integral.surface_integral, integral.surface_integral.integrand_hessian_u, Uu, p, Vv, domain)
end

function value!(val, integral::VariationalIntegral, Uu, p, domain)
  integrate!(val, integral.volume_integral, integral.volume_integral.integrand, Uu, p, domain)
  integrate!(val, integral.surface_integral, integral.surface_integral.integrand, Uu, p, domain)
end

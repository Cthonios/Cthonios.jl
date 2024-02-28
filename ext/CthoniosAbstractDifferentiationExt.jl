module CthoniosAbstractDifferentiationExt

import AbstractDifferentiation as AD
using Cthonios

function Cthonios.energy_gradient_u(domain::Cthonios.QuasiStaticDomain, Uu::V, backend) where V <: AbstractVector
  R = AD.gradient(backend, x -> Cthonios.energy(domain, x, domain.coords), Uu)[1]
  return R
end

function Cthonios.energy_gradient_x(domain::Cthonios.QuasiStaticDomain, Uu::V, backend) where V <: AbstractVector
  R = AD.gradient(backend, x -> Cthonios.energy(domain, Uu, x), domain.coords)[1]
  return R
end

hvp_u_func(domain::Cthonios.QuasiStaticDomain, Uu::V, backend) where V <: AbstractVector =
AD.pushforward_function(backend, x -> Cthonios.energy_gradient_u(domain, x, backend), Uu)

function Cthonios.energy_hvp_u(domain::Cthonios.QuasiStaticDomain, Uu::V, v::V, backend) where V <: AbstractVector
  func = hvp_u_func(domain, Uu, backend)
  return func(v)[1]
end

end # module
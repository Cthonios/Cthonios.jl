module CthoniosEnzymeExt

using Cthonios
using Enzyme
using EnzymeCore
using LinearAlgebra

function Cthonios.residual!(
  W::V1, dW::V1,
  domain::Cthonios.QuasiStaticDomain, 
  Uu::V1, dUu::V1, 
  U::V2, dU::V2,
  ::EnzymeCore.ReverseMode
) where {V1, V2}
  # autodiff_deferred(Reverse, energy!, Duplicated(W, dW), domain, Duplicated(Uu, dUu))

  autodiff_deferred(
    Reverse, Cthonios.energy!,
    Duplicated(W, dW),
    domain,
    Duplicated(Uu, dUu),
    Duplicated(U, dU)
    # U
  )
end

function residual_dot_v!(
  r_dot_v::V,
  W::V, dW::V,
  domain::Cthonios.QuasiStaticDomain, 
  Uu::V, dUu::V,
  v::V,
  mode::EnzymeCore.ReverseMode
) where V <: AbstractVector

  Cthonios.residual!(W, dW, domain, Uu, dUu, mode)
  r_dot_v[1] = dot(dUu, v)
  nothing
end

# function Cthonios.residual(domain::Cthonios.QuasiStaticDomain, Uu::V, mode::EnzymeCore.ReverseMode) where V <: AbstractVector
#   dUu = zeros(eltype(Uu), length(Uu))
#   W = [0.0]
#   dW = [1.0]
#   Cthonios.residual!(W, dW, domain, Uu, dUu, mode)
#   return dUu
# end

function Cthonios.residual(
  domain::Cthonios.QuasiStaticDomain, 
  Uu::V1, U::V2, 
  mode::EnzymeCore.ReverseMode
) where {V1, V2}
  dU = similar(U)
  dU .= 0.0
  dUu = zeros(eltype(Uu), length(Uu))
  W = [0.0]
  dW = [1.0]
  Cthonios.residual!(W, dW, domain, Uu, dUu, U, dU, mode)
  return dUu
end

function Cthonios.hvp_energy_u(
  domain::Cthonios.QuasiStaticDomain, 
  Uu::V1, U::V2
) where {V1, V2}

  # primals
  W = [0.0]

  # Reverse
  bW = [1.0]
  bU = similar(U)
  bU .= 0.0
  bUu = zeros(eltype(Uu), length(Uu))
  
  # Forward
  dW = [0.0]
  dUu = zeros(eltype(Uu), length(Uu))
  dUu[1] = 1.0
  dU = similar(U)
  dU .= 0.0

  # Forward over Reverse
  dbW = [0.0]
  dbUu = zeros(eltype(Uu), length(Uu))
  dbU = similar(U)
  U .= 0.0

  autodiff(
    Forward,
    (x, y, z) -> autodiff_deferred(Reverse, Cthonios.energy!, x, domain, y, z),
    Duplicated(Duplicated(W, bW), Duplicated(dW, dbW)),
    Duplicated(Duplicated(Uu, bUu), Duplicated(dUu, dbUu)),
    Duplicated(Duplicated(U, bU), Duplicated(dU, dbU)) 
  )
end

end # module
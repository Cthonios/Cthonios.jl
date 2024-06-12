module CthoniosEnzymeExt

using Cthonios
using Enzyme
using FiniteElementContainers

# helpers
function setup_dcache(::ReverseMode, cache::Cthonios.DomainCache)
  dcache = similar(cache)
  dcache.X .= zero(eltype(dcache.X))
  dcache.U .= zero(eltype(dcache.U))
  dcache.U_bc .= zero(eltype(dcache.U_bc))
  dcache.props .= zero(eltype(dcache.props))
  dcache.state_old .= zero(eltype(dcache.state_old))
  dcache.state_new .= zero(eltype(dcache.state_new))
  dcache.solver_cache.Π .= one(eltype(dcache.solver_cache.Π))
  dcache.solver_cache.Πs .= zero(eltype(dcache.solver_cache.Πs))
  dcache.solver_cache.assembler.residuals .= 
    zero(eltype(dcache.solver_cache.assembler.residuals))
  dcache.solver_cache.assembler.stiffnesses .= 
    zero(eltype(dcache.solver_cache.assembler.stiffnesses))
  return dcache
end

function setup_dUu(::ReverseMode, Uu)
  dUu = similar(Uu)
  dUu .= zero(eltype(dUu))
  return dUu
end

# Cthonios methods
function Cthonios.gradient!(dcache, ::ReverseMode, static, cache)
  autodiff(
    Reverse, Cthonios.internal_energy!,
    Const(static),
    Duplicated(cache, dcache)
  )
  return nothing
end

function Cthonios.gradient(::ReverseMode, obj::Objective)
  static, cache = obj.domain.static, obj.domain.cache
  dcache = setup_dcache(Reverse, cache)
  Cthonios.gradient!(dcache, Reverse, static, cache)
  return dcache
end

function Cthonios.internal_force!(dcache, ::ReverseMode, static, cache)
  autodiff(
    Reverse, Cthonios.internal_energy!, 
    Const(static),
    Duplicated(cache, dcache)
  )
  return nothing
end

function Cthonios.internal_force(::ReverseMode, domain::Domain)
  static, cache = domain.static, domain.cache
  dcache = setup_dcache(Reverse, cache)
  Cthonios.internal_force!(dcache, Reverse, static, cache)
  return dcache.U
end

function Cthonios.residual!(dcache, dUu, ::ReverseMode, static, cache, Uu)
  autodiff(
    Reverse, Cthonios.internal_energy!, 
    Const(static), 
    Duplicated(cache, dcache),
    Duplicated(Uu, dUu)
  )
  return nothing
end

function Cthonios.residual(::ReverseMode, domain::Domain, Uu)
  static, cache = domain.static, domain.cache
  dcache = setup_dcache(Reverse, cache)
  dUu = setup_dUu(Reverse, Uu)
  Cthonios.residual!(dcache, dUu, Reverse, static, cache, Uu)
  return dUu
end

function Cthonios.residual(::ReverseMode, obj::Cthonios.Objective, Uu)
  static, cache = obj.domain.static, obj.domain.cache
  dcache = setup_dcache(Reverse, cache)
  dUu = setup_dUu(Reverse, Uu)
  autodiff(
    Reverse,
    obj.value,
    Const(static),
    Duplicated(cache, dcache),
    Duplicated(Uu, dUu)
  )
  dUu
end

end # module
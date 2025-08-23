module CthoniosAdaptExt

using Adapt
using Cthonios
using FiniteElementContainers

function Adapt.adapt(to, cache::Cthonios.QuadratureLevelObjectiveCache)
  return Cthonios.QuadratureLevelObjectiveCache(
    cache.objective,
    adapt(to, cache.sim_cache),
    cache.timer
  )
end

function Adapt.adapt(to, cache::Cthonios.SingleDomainSimulationCache)
  return Cthonios.SingleDomainSimulationCache(
    adapt(to, cache.assembler),
    adapt(to, cache.parameters),
    cache.post_processor,
    cache.time
  )
end

end # module

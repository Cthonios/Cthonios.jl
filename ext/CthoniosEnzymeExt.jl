module CthoniosEnzymeExt

using Cthonios
using Enzyme

function Cthonios.energy_gradient(solver, domain, backend)

  # unpack primal cache and make a copy for the tangent
  primal_cache = domain.domain_cache
  btangent_cache = similar(domain.domain_cache)

  # unpack specific arrays
  Uu = solver.Uu
  dUu = similar(Uu)
  dUu .= 0.0
  
  # seed grads
  btangent_cache.X .= 0.0
  btangent_cache.U .= 0.0
  btangent_cache.props .= 0.0
  btangent_cache.state_old .= 0.0
  btangent_cache.state_new .= 0.0
  btangent_cache.Π .= 1.0
  btangent_cache.Πs .= 1.0
  btangent_cache.f .= 0.0
  btangent_cache.V .= 0.0

  autodiff(
    Reverse, Cthonios.energy_new!,
    Const(solver),
    Const(domain),
    Duplicated(Uu, dUu),
    Duplicated(primal_cache, btangent_cache),
    Const(backend)
  )

  btangent_cache
end

function grad!(solver, domain, Uu, cache, backend)
  autodiff_deferred(
    Reverse, Cthonios.energy_new!,
    solver, domain,
    Uu, cache,
    backend
  )
  return nothing
end

# still not working
function Cthonios.energy_hvp(solver, domain, backend)
  # unpack primal cache and make a copy for the tangent
  primal_cache = domain.domain_cache
  btangent_cache = similar(domain.domain_cache)
  dtangent_cache = similar(domain.domain_cache)
  dbtangent_cache = similar(domain.domain_cache)

  # unpack specific arrays
  Uu = solver.Uu
  bUu = similar(Uu)
  dUu = similar(Uu)
  dbUu = similar(Uu)

  # seed backward grads
  btangent_cache.X .= 0.0
  btangent_cache.U .= 0.0
  btangent_cache.props .= 0.0
  btangent_cache.state_old .= 0.0
  btangent_cache.state_new .= 0.0
  btangent_cache.Π .= 1.0
  btangent_cache.Πs .= 1.0
  btangent_cache.f .= 0.0
  btangent_cache.V .= 0.0
  bUu .= 0.0

  # seed forward grads
  dtangent_cache.X .= 0.0
  dtangent_cache.U .= 0.0
  dtangent_cache.props .= 1.0
  dtangent_cache.state_old .= 1.0
  dtangent_cache.state_new .= 1.0
  dtangent_cache.Π .= 0.0
  dtangent_cache.Πs .= 0.0
  dtangent_cache.f .= 1.0
  dtangent_cache.V .= 1.0
  dUu .= 1.0

  # seed forward over reverse grads
  dbtangent_cache.X .= 0.0
  dbtangent_cache.U .= 0.0
  dbtangent_cache.props .= 0.0
  dbtangent_cache.state_old .= 0.0
  dbtangent_cache.state_new .= 0.0
  dbtangent_cache.Π .= 0.0
  dbtangent_cache.Πs .= 0.0
  dbtangent_cache.f .= 0.0
  dbtangent_cache.V .= 0.0
  dbUu .= 0.0

  # grad!(solver, domain, Duplicated(Uu, dUu), Duplicated(primal_cache, btangent_cache), backend)

  autodiff(
    Forward, grad!,
    Const(solver),
    Const(domain),
    Duplicated(
      Duplicated(Uu, bUu), 
      Duplicated(dUu, dbUu)
    ),
    Duplicated(
      Duplicated(primal_cache, btangent_cache), 
      Duplicated(dtangent_cache, dbtangent_cache),
    ), 
    Const(backend)
  )


  # btangent_cache
end

end # module
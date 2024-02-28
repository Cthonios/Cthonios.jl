module CthoniosAdaptExt

using Adapt
using Cthonios

function Adapt.adapt_structure(to, section::Cthonios.TotalLagrangeSection)
  fspace = Adapt.adapt_structure(to, section.fspace)
  formulation = Adapt.adapt_structure(to, section.formulation)
  model = Adapt.adapt_structure(to, section.model)
  return Cthonios.TotalLagrangeSection(fspace, formulation, model)
end

function Adapt.adapt_structure(to, cache::Cthonios.QuasiStaticDomainCache)
  X = Adapt.adapt_structure(to, cache.X)
  U = Adapt.adapt_structure(to, cache.U)
  state = Adapt.adapt_structure(to, cache.state)
  props = Adapt.adapt_structure(to, cache.props)
  Π = Adapt.adapt_structure(to, cache.Π)
  Πs = Adapt.adapt_structure(to, cache.Πs)
  V = Adapt.adapt_structure(to, cache.V)
  return Cthonios.QuasiStaticDomainCache(X, U, state, props, Π, Πs, V)
end

function Adapt.adapt_structure(to, domain::Cthonios.QuasiStaticDomain)
  dof = Adapt.adapt_structure(to, domain.dof)
  # funcs = Adapt.adapt_structure(to, domain.funcs)
  funcs = domain.funcs # how to put this on a GPU?
  bc_nodes = Adapt.adapt_structure(to, domain.bc_nodes)
  bc_dofs = Adapt.adapt_structure(to, domain.bc_dofs)
  bc_func_ids = Adapt.adapt_structure(to, domain.bc_func_ids)
  sections = Adapt.adapt_structure(to, domain.sections)
  time = Adapt.adapt_structure(to, domain.time)
  domain_cache = Adapt.adapt_structure(to, domain.domain_cache)
  return Cthonios.QuasiStaticDomain(dof, funcs, bc_nodes, bc_dofs, bc_func_ids, sections, time, domain_cache)
end

function Adapt.adapt_structure(to, solver::Cthonios.DirectLinearSolver)
  settings = Adapt.adapt_structure(to, solver.settings)
  assembler = Adapt.adapt_structure(to, solver.assembler)
  factorization = Adapt.adapt_structure(to, solver.factorization)
  return Cthonios.DirectLinearSolver(settings, assembler, factorization)
end

function Adapt.adapt_structure(to, solver::Cthonios.NewtonSolver)
  settings = Adapt.adapt_structure(to, solver.settings)
  linear_solver = Adapt.adapt_structure(to, solver.linear_solver)
  Uu = Adapt.adapt_structure(to, solver.Uu)
  ΔUu = Adapt.adapt_structure(to, solver.ΔUu)
  return Cthonios.NewtonSolver(settings, linear_solver, Uu, ΔUu)
end

function Adapt.adapt_structure(to, solver::Cthonios.TrustRegionSolver)
  settings = Adapt.adapt_structure(to, solver.settings)
  linear_solver = Adapt.adapt_structure(to, solver.linear_solver)
  Uu = Adapt.adapt_structure(to, solver.Uu)
  ΔUu = Adapt.adapt_structure(to, solver.ΔUu)
  Hv = Adapt.adapt_structure(to, solver.Hv)
  Hv_field = Adapt.adapt_structure(to, solver.Hv_field)
  return Cthonios.TrustRegionSolver(settings, linear_solver, Uu, ΔUu, Hv, Hv_field)
end

function Adapt.adapt_structure(to, prob::Cthonios.ForwardProblem)
  domain = Adapt.adapt_structure(to, prob.domain)
  solver = Adapt.adapt_structure(to, prob.solver)
  post_processor = Adapt.adapt_structure(to, prob.post_processor)
  return Cthonios.ForwardProblem(domain, solver, post_processor)
end

end # module
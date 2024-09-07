"""
$(TYPEDEF)
"""
abstract type AbstractObjective end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Objective{D, F1, F2, F3}
  domain::D
  value::F1
  gradient::F2
  hessian::F3
end

"""
$(TYPEDEF)
"""
abstract type AbstractObjectiveParameters end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct ObjectiveParameters{U1, T, U2, U3} <: AbstractObjectiveParameters
  X::U1
  t::T
  U::U2
  hvp_scratch::U3
end

function ObjectiveParameters(o::Objective, times)
  X = coordinates(o.domain.mesh)
  X = NodalField{size(X), Vector}(X)
  U = create_fields(o.domain)
  hvp_scratch = create_fields(o.domain)
  params = ObjectiveParameters(X, times, U, hvp_scratch)
  return params
end

function gradient!(g, o::Objective, Uu, p::ObjectiveParameters)
  g .= zero(eltype(g))
  update_fields!(p.U, o.domain, Uu)
  domain_iterator!(g, o.gradient, o.domain, p.U, p.X)
  return @views g[o.domain.dof.unknown_dofs]
end

function hessian!(asm::FiniteElementContainers.StaticAssembler, o::Objective, Uu, p::ObjectiveParameters)
  asm.stiffnesses .= zero(eltype(asm.stiffnesses))
  update_fields!(p.U, o.domain, Uu)
  domain_iterator!(asm, o.hessian, o.domain, p.U, p.X)
  return SparseArrays.sparse!(asm) |> Symmetric
end

function hvp!(Hv, o::Objective, Uu, p::ObjectiveParameters, Vv)
  Hv .= zero(eltype(Hv))
  update_fields!(p.U, o.domain, Uu)
  update_fields!(p.hvp_scratch, o.domain, Vv)
  domain_iterator!(Hv, o.hessian, o.domain, p.U, p.X, p.hvp_scratch)
  return @views Hv[o.domain.dof.unknown_dofs]
end

function objective!(val, o::Objective, Uu, p::ObjectiveParameters)
  val .= zero(eltype(val))
  update_fields!(p.U, o.domain, Uu)
  domain_iterator!(val, o.value, o.domain, p.U, p.X)
  return @views val[1]
end

function step!(p::ObjectiveParameters)
  step!(p.t)
  return nothing
end

function update_bcs!(p::ObjectiveParameters, o::Objective)
  update_bcs!(p.U, o.domain, p.X, p.t)
  return nothing
end

# exports
export Objective
export ObjectiveParameters

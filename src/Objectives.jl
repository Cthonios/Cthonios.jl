"""
$(TYPEDEF)
"""
abstract type AbstractObjective end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Objective{D, F1, F2, F3, T}
  domain::D
  value::F1
  gradient::F2
  hessian::F3
  timer::T
end

timer(o::Objective) = o.timer

"""
$(TYPEDEF)
"""
abstract type AbstractObjectiveParameters end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct ObjectiveParameters{U1, T, B, U2, U3} <: AbstractObjectiveParameters
  # design parameters
  X::U1
  t::T
  Ubc::B
  # scratch arrays
  U::U2
  hvp_scratch::U3
end

function ObjectiveParameters(o::Objective, times)
  X = coordinates(o.domain.mesh)
  X = NodalField{size(X), Vector}(X)
  U = create_fields(o.domain)
  # boundary conditions
  Ubc = Vector{eltype(X)}(undef, 0)
  update_dirichlet_vals!(Ubc, o.domain, X, times)
  # scratch arrays
  hvp_scratch = create_fields(o.domain)
  params = ObjectiveParameters(
    X, times, Ubc, 
    U, hvp_scratch
  )
  return params
end

function gradient!(g, o::Objective, Uu, p::ObjectiveParameters)
  g .= zero(eltype(g))
  update_field_unknowns!(p.U, o.domain, Uu)
  domain_iterator!(g, o.gradient, o.domain, Uu, p)
  return @views g[o.domain.dof.unknown_dofs]
end

# used in EnzymeExt
function gradient!(g, o::Objective, U, X::NodalField)
  @timeit timer(o) "Objective - gradient!" begin
    g .= zero(eltype(g))
    domain_iterator!(g, o.gradient, o.domain, U, X)
  end
  return nothing
end

function hessian!(asm::FiniteElementContainers.StaticAssembler, o::Objective, Uu, p::ObjectiveParameters)
  @timeit timer(o) "Objective - hessian!" begin
    asm.stiffnesses .= zero(eltype(asm.stiffnesses))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(asm, o.hessian, o.domain, Uu, p)
    return SparseArrays.sparse!(asm) |> Symmetric
  end
end

function hvp!(Hv, o::Objective, Uu, p::ObjectiveParameters, Vv)
  @timeit timer(o) "Objective - hvp!" begin
    Hv .= zero(eltype(Hv))
    update_field_unknowns!(p.U, o.domain, Uu)
    update_field_unknowns!(p.hvp_scratch, o.domain, Vv)
    domain_iterator!(Hv, o.hessian, o.domain, Uu, p, Vv)
    return @views Hv[o.domain.dof.unknown_dofs]
  end
end

function objective(o::Objective, U::NodalField, Ubc, X::NodalField)
  val = zeros(1)
  domain_iterator!(val, o.value, o.domain, U, Ubc, X)
  return val
end

function objective!(val, o::Objective, Uu, p::ObjectiveParameters)
  @timeit timer(o) "Objective - objective!" begin
    val .= zero(eltype(val))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(val, o.value, o.domain, Uu, p)
  end
  return @views val[1]
end

function step!(p::ObjectiveParameters)
  step!(p.t)
  return nothing
end

function update_dirichlet_vals!(p::ObjectiveParameters, o::Objective)
  update_dirichlet_vals!(p.Ubc, o.domain, p.X, p.t)
  return nothing
end

# exports
export Objective
export ObjectiveParameters

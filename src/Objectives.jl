"""
$(TYPEDEF)
Abstract base objective type.
"""
abstract type AbstractObjective end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Objective type for evaluating objective functions.
The functions correspond to quadrature level evaluations
of the objective function, it's gradient, and it's hessian
respectively.
```domain``` - A domain object 
```value``` - A function for the quadrature level objective function evaluation.
```gradient``` - A function for the quadrature level objective function gradient evaluation.
```hessian``` - A function for the quadrature level objective hessian evaluation.
```timer``` - A timer that's already setup.
"""
struct Objective{D, F1, F2, F3, F4, F5, F6, T}
  domain::D
  value::F1
  gradient::F2
  hessian::F3
  neumann_value::F4
  neumann_gradient::F5
  neumann_hessian::F6
  timer::T
end

function Objective(domain::Domain, value_func::F1, neumann_value_func::F2, timer::TimerOutput) where {F1 <: Function, F2 <: Function}
  grad_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> value_func(phys, cell, z, x, state, props, t), u)
  hess_func(phys, cell, u, x, state, props, t) = ForwardDiff.hessian(z -> value_func(phys, cell, z, x, state, props, t), u)
  neumann_grad_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_value_func(phys, cell, z, x, t, bc, val, r, q, f), u)
  neumann_hess_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.hessian(z -> neumann_value_func(phys, cell, z, x, t, bc, val, r, q, f), u)
  return Objective(
    domain, 
    value_func, grad_func, hess_func, 
    neumann_value_func, neumann_grad_func, neumann_hess_func,
    timer
  )
end

function Objective(domain::Domain, ::typeof(energy))
  return Objective(
    domain,
    energy, gradient, hessian,
    neumann_energy, neumann_gradient, neumann_hessian,
    TimerOutput()
  )
end

function Objective(inputs::Dict{Symbol, Any}, domain, timer)
  value_func = eval(Meta.parse(inputs[:value]))
  grad_func = eval(Meta.parse(inputs[:gradient]))
  hess_func = eval(Meta.parse(inputs[:hessian]))
  neumann_value_func = eval(Meta.parse(inputs[Symbol("neumann value")]))
  neumann_grad_func = eval(Meta.parse(inputs[Symbol("neumann gradient")]))
  neumann_hess_func = eval(Meta.parse(inputs[Symbol("neumann hessian")]))
  return Objective(
    domain, 
    value_func, grad_func, hess_func, 
    neumann_value_func, neumann_grad_func, neumann_hess_func, 
    timer
  )
end

function grad_u end
function grad_p end
function val_and_grad_u end
function val_and_grad_p end
function hvp_u end
function hvp_p end

timer(o::Objective) = o.timer

"""
$(TYPEDEF)
Abstract base type for objective function parameters.
"""
abstract type AbstractObjectiveParameters end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Type for objective function parameters for design parameters
such as coordinates, time, bc values, properties
state variables, and some scratch arrays.
"""
struct ObjectiveParameters{U1, T, B, N, S, P, U2, U3, Q} <: AbstractObjectiveParameters
  # design parameters
  X::U1
  t::T
  Ubc::B # dirichlet bc design parameters
  nbc::N # neumann bc design parameters
  state_old::S
  state_new::S
  props::P
  # scratch arrays
  U::U2
  hvp_scratch::U3
  q_vals_scratch::Q
end

"""
$(TYPEDSIGNATURES)
Constructor for a ```ObjectiveParameters``` type.
```o``` - Objective function object.
```times``` - Times object.
"""
function ObjectiveParameters(o::Objective, times)
  X = copy(o.domain.coords)
  U = create_fields(o.domain)
  # boundary conditions
  Ubc = Vector{eltype(X)}(undef, 0)
  nbc = Vector{SVector{size(X, 1), eltype(X)}}(undef, 0)
  update_dirichlet_vals!(Ubc, o.domain, X, times)
  update_neumann_vals!(nbc, o.domain, X, times)
  # properties
  props = map(sec -> sec.props, o.domain.sections)
  props = ComponentArray(props)
  # scratch arrays
  hvp_scratch = create_fields(o.domain)
  # need a scratch array for calculating q values on gpus
  # TODO move somewhere else
  q_vals_scratch = Dict{Symbol, Any}()
  state_old = Dict{Symbol, Any}()
  state_new = Dict{Symbol, Any}()
  for (name, sec) in pairs(o.domain.sections)
    NQ = FiniteElementContainers.num_q_points(sec.fspace)
    NE = FiniteElementContainers.num_elements(sec.fspace)
    q_vals_scratch[name] = zeros(eltype(X), NQ, NE)
    state_old[name] = repeat(init_state(sec.physics), outer=(1, NQ, NE))
    state_new[name] = repeat(init_state(sec.physics), outer=(1, NQ, NE))
  end
  q_vals_scratch = ComponentArray(q_vals_scratch)
  state_old = ComponentArray(state_old)
  state_new = ComponentArray(state_new)
  params = ObjectiveParameters(
    X, times, Ubc, nbc, state_old, state_new, props,
    U, hvp_scratch, q_vals_scratch
  )
  return params
end

function gradient(o::Objective, Uu, p)
  @timeit timer(o) "Objective - gradient" begin
    g = similar(p.hvp_scratch)
    g .= zero(eltype(g))
    gradient!(g, o, Uu, p)
    return g[o.domain.dof.unknown_dofs]
  end
end

"""
$(TYPEDSIGNATURES)
"""
function gradient!(g, o::Objective, Uu, p::ObjectiveParameters)
  @timeit timer(o) "Objectives - gradient!" begin
    g .= zero(eltype(g))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(g, o.gradient, o.domain, Uu, p)
    surface_iterator!(g, o.neumann_gradient, o.domain, Uu, p)
    return @views g[o.domain.dof.unknown_dofs]
  end
end

function gradient_for_ad!(g, o::Objective, Uu, p::ObjectiveParameters)
  # @timeit timer(o) "Objectives - gradient!" begin
    g .= zero(eltype(g))

    # TTRYING BELOW OUT
    update_neumann_vals!(p, o)
    # TRYING ABOVE out

    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(g, o.gradient, o.domain, Uu, p)
    # surface_iterator!(g, o.neumann_gradient, o.domain, Uu, p)
    # return @views g[o.domain.dof.unknown_dofs]
    return nothing
  # end
end

"""
$(TYPEDSIGNATURES)
"""
function hessian!(asm::FiniteElementContainers.Assembler, o::Objective, Uu, p::ObjectiveParameters)
  @timeit timer(o) "Objective - hessian!" begin
    asm.stiffnesses .= zero(eltype(asm.stiffnesses))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(asm, o.hessian, o.domain, Uu, p)
    # surface_iterator!(asm, o.neumann_hessian, o.domain, Uu, p) # Don't even need this except for maybe some cases?
    return SparseArrays.sparse!(asm) |> Symmetric
  end
end

function hvp(o::Objective, Uu, p, Vv)
  @timeit timer(o) "Objective - hvp" begin
    Hv = similar(p.hvp_scratch)
    Hv .= zero(eltype(Hv))
    hvp!(Hv, o, Uu, p, Vv)
    return @views Hv[o.domain.dof.unknown_dofs]
  end
end

"""
$(TYPEDSIGNATURES)
"""
function hvp!(Hv::AbstractVector, o::Objective, Uu, p::ObjectiveParameters, Vv)
  @timeit timer(o) "Objective - hvp!" begin
    Hv .= zero(eltype(Hv))
    hvp!(p.hvp_scratch, o, Uu, p, Vv)
    @views Hv .= p.hvp_scratch[o.domain.dof.unknown_dofs]
  end
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function hvp!(Hv::NodalField, o::Objective, Uu, p::ObjectiveParameters, Vv)
  @timeit timer(o) "Objective - hvp!" begin
    Hv .= zero(eltype(Hv))
    update_field_unknowns!(p.U, o.domain, Uu)
    update_field_unknowns!(p.hvp_scratch, o.domain, Vv)
    domain_iterator!(Hv, o.hessian, o.domain, Uu, p, Vv)
    # TODO do we need a surface iterator for more general neumann bcs?
    return @views Hv[o.domain.dof.unknown_dofs]
  end
end

# """
# $(TYPEDSIGNATURES)
# """
# function objective(o::Objective, Uu, p)
#   # return domain_iterator(o.value, o.domain, Uu, p)
#   return domain_iterator(energy, o.domain, Uu, p)
# end

# function grad_u(o::Objective, Uu, p, backend)
#   # func = x -> objective(o, x, p)
#   # DifferentiationInterface.gradient(func, backend, Uu)
#   # DifferentiationInterface.gradient(x -> objective(o, x, p), backend, Uu)
#   func = x -> domain_iterator(energy, o.domain, x, p)
#   DifferentiationInterface.gradient(func, backend, Uu)
# end

function mass_matrix!(asm, o::Objective, Uu, p::ObjectiveParameters)
  @timeit timer(o) "Objectives - mass_matrix!" begin
    # asm.masses .= zero(eltype(asm.masses))
    asm.stiffnesses .= zero(eltype(asm.stiffnesses))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(asm, mass_matrix, o.domain, Uu, p)
    return SparseArrays.sparse!(asm) |> Symmetric
  end
end

function objective(o::Objective, Uu, p)
  @timeit timer(o) "Objective - objective" begin
    val = zeros(1)
    objective!(val, o, Uu, p)
    return val[1]
  end
end

"""
$(TYPEDSIGNATURES)
"""
function objective!(val, o::Objective, Uu, p::ObjectiveParameters)
  @timeit timer(o) "Objective - objective!" begin
    val .= zero(eltype(val))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(val, o.value, o.domain, Uu, p)
    surface_iterator!(val, o.neumann_value, o.domain, Uu, p)
    return @views val[1]
  end
end

"""
$(TYPEDSIGNATURES)
"""
function step!(p::ObjectiveParameters)
  step!(p.t)
  return nothing
end

function step_new!(p::ObjectiveParameters, o::Objective)
  step!(p.t)
  update_dirichlet_vals!(p, o)
  update_neumann_vals!(p, o)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_dirichlet_vals!(p::ObjectiveParameters, o::Objective)
  update_dirichlet_vals!(p.Ubc, o.domain, p.X, p.t)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_neumann_vals!(p::ObjectiveParameters, o::Objective)
  update_neumann_vals!(p.nbc, o.domain, p.X, p.t)
  return nothing
end

# exports
export Objective
export ObjectiveParameters

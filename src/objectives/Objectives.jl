"""
$(TYPEDEF)
Abstract base objective type.
"""
abstract type AbstractObjective end
timer(o::AbstractObjective) = o.timer
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
# struct Objective{D, F1, F2, F3, F4, F5, F6, T} <: AbstractObjective
#   domain::D
#   value::F1
#   gradient::F2
#   hessian::F3
#   neumann_value::F4
#   neumann_gradient::F5
#   neumann_hessian::F6
#   timer::T
# end

# function Objective(domain::Domain, value_func::F1, neumann_value_func::F2, timer::TimerOutput) where {F1 <: Function, F2 <: Function}
#   grad_func(phys, cell, u, x, state, props, t) = ForwardDiff.gradient(z -> value_func(phys, cell, z, x, state, props, t), u)
#   hess_func(phys, cell, u, x, state, props, t) = ForwardDiff.hessian(z -> value_func(phys, cell, z, x, state, props, t), u)
#   neumann_grad_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.gradient(z -> neumann_value_func(phys, cell, z, x, t, bc, val, r, q, f), u)
#   neumann_hess_func(phys, cell, u, x, t, bc, val, r, q, f) = ForwardDiff.hessian(z -> neumann_value_func(phys, cell, z, x, t, bc, val, r, q, f), u)
#   return Objective(
#     domain, 
#     value_func, grad_func, hess_func, 
#     neumann_value_func, neumann_grad_func, neumann_hess_func,
#     timer
#   )
# end

# function Objective(domain::Domain, ::typeof(energy))
#   return Objective(
#     domain,
#     energy, gradient, hessian,
#     neumann_energy, neumann_gradient, neumann_hessian,
#     TimerOutput()
#   )
# end

# function Objective(inputs::Dict{Symbol, Any}, domain, timer)
#   value_func = eval(Meta.parse(inputs[:value]))

#   if typeof(value_func) == typeof(energy)
#     return Objective(domain, value_func)
#   else
#     throw(ErrorException("Unsupported objective function"))
#   end
# end

# function grad_u end
# function grad_p end
# function val_and_grad_u end
# function val_and_grad_p end
# function hvp_u end
# function hvp_p end

# timer(o::Objective) = o.timer

"""
$(TYPEDEF)
Abstract base type for objective function parameters.
"""
abstract type AbstractObjectiveParameters end


function gradient(o::AbstractObjective, Uu, p)
  @timeit timer(o) "Objective - gradient" begin
    g = p.grad_scratch
    g .= zero(eltype(g))
    gradient!(g, o, Uu, p)
    # update vec now
    g_vec = p.grad_vec_scratch
    g_vec .= @views g.vals[o.domain.dof.unknown_dofs]
    return g_vec
  end
end

"""
$(TYPEDSIGNATURES)
"""
function gradient!(g, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  @timeit timer(o) "Objectives - gradient!" begin
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(g, o.gradient, o.domain, Uu, p)
    surface_iterator!(g, o.neumann_gradient, o.domain, Uu, p)
  end
end

function gradient_for_ad!(g, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  # @timeit timer(o) "Objectives - gradient!" begin
    g .= zero(eltype(g))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(g, o.gradient, o.domain, Uu, p)
    surface_iterator!(g, o.neumann_gradient, o.domain, Uu, p)
    # return @views g[o.domain.dof.unknown_dofs]
    return nothing
  # end
end

"""
$(TYPEDSIGNATURES)
"""
function hessian!(asm::FiniteElementContainers.Assembler, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  @timeit timer(o) "Objective - hessian!" begin
    asm.stiffnesses .= zero(eltype(asm.stiffnesses))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(asm, o.hessian, o.domain, Uu, p)
    # surface_iterator!(asm, o.neumann_hessian, o.domain, Uu, p) # Don't even need this except for maybe some cases?
    return SparseArrays.sparse!(asm) |> Symmetric
  end
end

# TODO need to make below have pre-allocated arrays
# currently breaks trust region solver though
# probably need to modularize that guy.
function hvp(o::AbstractObjective, Uu, p, Vv)
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
function hvp!(Hv::AbstractVector, o::AbstractObjective, Uu, p::AbstractObjectiveParameters, Vv)
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
function hvp!(Hv::NodalField, o::AbstractObjective, Uu, p::AbstractObjectiveParameters, Vv)
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

function mass_matrix!(asm, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  @timeit timer(o) "Objectives - mass_matrix!" begin
    # asm.masses .= zero(eltype(asm.masses))
    asm.stiffnesses .= zero(eltype(asm.stiffnesses))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(asm, mass_matrix, o.domain, Uu, p)
    return SparseArrays.sparse!(asm) |> Symmetric
  end
end

function objective(o::AbstractObjective, Uu, p)
  @timeit timer(o) "Objective - objective" begin
    val = zeros(1)
    objective!(val, o, Uu, p)
    return val[1]
  end
end

"""
$(TYPEDSIGNATURES)
"""
function objective!(val, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  @timeit timer(o) "Objective - objective!" begin
    val .= zero(eltype(val))
    update_field_unknowns!(p.U, o.domain, Uu)
    domain_iterator!(val, o.value, o.domain, Uu, p)
    surface_iterator!(val, o.neumann_value, o.domain, Uu, p)
    return @views val[1]
  end
end

# implementations
include("ObjectiveParameters.jl")
include("UnconstrainedObjective.jl")

# exports
# export Objective
export ObjectiveParameters
export UnconstrainedObjective

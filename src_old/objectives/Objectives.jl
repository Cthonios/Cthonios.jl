"""
$(TYPEDEF)
Abstract base objective type.
"""
abstract type AbstractObjective end
timer(o::AbstractObjective) = o.timer

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
    update_field_unknowns!(current_solution(p), o.domain, Uu)
    gradient_u!(g, o.integral, Uu, p, o.domain)
  end
end

function gradient_x(o::AbstractObjective, Uu, p)
  @timeit timer(o) "Objective - gradient_x" begin
    dX = similar(p.X)
    dX .= zero(eltype(dX))
    gradient_x!(dX, o, Uu, p)
    return dX
  end
end

function gradient_x!(g, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  @timeit timer(o) "Objectives - gradient_x!" begin
    update_field_unknowns!(current_solution(p), o.domain, Uu)
    gradient_x!(g, o.integral, Uu, p, o.domain)
  end
end

function gradient_for_ad!(g, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  # @timeit timer(o) "Objectives - gradient!" begin
    g .= zero(eltype(g))
    update_field_unknowns!(current_solution(p), o.domain, Uu)
    gradient_u!(g, o.integral, Uu, p, o.domain)
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
    update_field_unknowns!(current_solution(p), o.domain, Uu)
    hessian_u!(asm, o.integral, Uu, p, o.domain)
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
    update_field_unknowns!(current_solution(p), o.domain, Uu)
    update_field_unknowns!(p.hvp_scratch, o.domain, Vv)
    hvp_u!(Hv, o.integral, Uu, p, Vv, o.domain)
    # TODO do we need a surface iterator for more general neumann bcs?
    return @views Hv[o.domain.dof.unknown_dofs]
  end
end

"""
$(TYPEDSIGNATURES)
Special case for use in dynamics problems
"""
function mass_matrix!(asm, o::AbstractObjective, Uu, p::AbstractObjectiveParameters)
  @timeit timer(o) "Objectives - mass_matrix!" begin
    # asm.masses .= zero(eltype(asm.masses))
    asm.stiffnesses .= zero(eltype(asm.stiffnesses))
    update_field_unknowns!(current_solution(p), o.domain, Uu)
    integrate!(asm, o.integral.volume_integral, mass_matrix, Uu, p, o.domain)
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
    update_field_unknowns!(current_solution(p), o.domain, Uu)
    value!(val, o.integral, Uu, p, o.domain)
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

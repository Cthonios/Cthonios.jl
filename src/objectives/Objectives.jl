abstract type AbstractObjective{
    A  <: FiniteElementContainers.AbstractAssembler,
    F1 <: Function,
    T  <: TimerOutput
} end

function FiniteElementContainers.create_field(o::AbstractObjective)
    return FiniteElementContainers.create_field(o.assembler, H1Field)
end

function FiniteElementContainers.create_unknowns(o::AbstractObjective)
    return FiniteElementContainers.create_unknowns(o.assembler, H1Field)
end

function gradient(o::AbstractObjective, Uu, p)
    @timeit o.timer "Objective - gradient!" begin
        assemble_vector!(o.assembler, Uu, p, H1Field, o.gradient_u)
    end
    return residual(o.assembler)
end

function hessian(o::AbstractObjective, Uu, p)
    @timeit o.timer "Objective - gradient!" begin
        assemble_matrix!(o.assembler, Uu, p, H1Field, o.hessian_u)
    end
    return stiffness(o.assembler) |> Symmetric # needed for cholesky
end

function hvp(o::AbstractObjective, Uu, p, Vu)
    @timeit o.timer "Objective - gradient!" begin
        # this first one allocates for some reason
        # assemble_matrix_action!(o.assembler, Uu, p, Vu, H1Field, o.hessian_u)
        # this one does not
        assemble!(o.assembler, Uu, p, Vu, Val{:stiffness_action}(), H1Field)
    end
    return FiniteElementContainers.hvp(o.assembler)
end

function value(o::AbstractObjective, Uu, p)
    @timeit o.timer "Objective - value!" begin
        assemble_scalar!(o.assembler, Uu, p, H1Field, o.value)
    end
    val = mapreduce(x -> sum(x), sum, o.assembler.scalar_quadarature_storage)
    return val
end

include("UnconstrainedObjective.jl")

"""
"""
abstract type AbstractObjective{
    RT <: Number,
    RV <: AbstractArray{RT, 1},
    A # Assembler type
} end

function FiniteElementContainers.create_field(o::AbstractObjective)
    return FiniteElementContainers.create_field(o.assembler)
end

function FiniteElementContainers.create_unknowns(o::AbstractObjective)
    return FiniteElementContainers.create_unknowns(o.assembler)
end

function assembler(o::AbstractObjective)
    return o.assembler
end

abstract type AbstractSolidMechanicsObjective{RT, RV, A} <: AbstractObjective{RT, RV, A} end

function _setup_solid_mechanics_assembler(
    mesh,
    q_degree = 2,
    use_condensed = false,
    use_inplace_methods = true
)
    fspace = FunctionSpace(mesh, H1Field, Lagrange)
    u = VectorFunction(fspace, "displ")
    dof = DofManager(u; use_condensed = use_condensed)
    assembler = SparseMatrixAssembler(dof; use_inplace_methods = use_inplace_methods)
    return assembler
end

# include("ContactObjective.jl")
# include("ConstrainedObjective.jl")
# include("DesignObjective.jl")
include("ExplicitDynamicsObjective.jl")
include("ImplicitDynamicsObjective.jl")
include("QuasiStaticObjective.jl")

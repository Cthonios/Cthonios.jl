Base.@kwdef struct OutputSettings
    # nodal variables
    acceleration::Bool       = false
    displacement::Bool       = true
    temperature::Bool        = false
    velocity::Bool           = false
    # element/quadrature variables
    cauchy_stress::Bool      = false
    # global variables
    external_energy::Bool    = false
    internal_energy::Bool    = false
    kinetic_energy::Bool     = false
    # output frequency
    output_exodus_every::Int = 1
end

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

function default_output_settings(::AbstractObjective)
    return OutputSettings()
end

# TODO need a lot of work of element variables
function initialize_postprocessor(mesh, output_file, objective, output_settings)
    fspace = FEC.function_space(assembler(objective))
    pp = PostProcessor(mesh, output_file)

    try
        if output_settings.acceleration
            a = VectorFunction(fspace, "acceleration")
            FEC.add_function!(pp, a)
        end

        if output_settings.displacement
            u = VectorFunction(fspace, "displ")
            FEC.add_function!(pp, u)
        end

        if output_settings.velocity
            v = VectorFunction(fspace, "velocity")
            FEC.add_function!(pp, v)
        end
    catch e
        close(pp)
    end

    FEC.finalize_setup!(pp)
    return pp
end

function run!(
    solver, objective::AbstractObjective, u, p,
    mesh, output_file, 
    output_settings = default_output_settings(objective);
    verbose::Bool = true
)
    pp = initialize_postprocessor(mesh, output_file, objective, output_settings)
    initialize!(objective, u, p)
    time_end = p.times.time_end
    n = 2
    try
        # post-process first step in try block
        postprocess!(pp, output_settings, n, objective, u, p)
        while FEC.current_time(p.times) < time_end - 1e3 * eps(time_end)
            step!(solver, objective, u, p; verbose = verbose)
            if n % output_settings.output_exodus_every == 0
                postprocess!(pp, output_settings, n, objective, u, p)
            end
            n = n + 1
        end
    finally
        close(pp)
    end
    return nothing
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
include("EigenObjective.jl")
include("ExplicitDynamicsObjective.jl")
include("ImplicitDynamicsObjective.jl")
include("QuasiStaticObjective.jl")

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

        if output_settings.temperature
            theta = ScalarFunction(fspace, "temperature")
            FEC.add_function!(pp, theta)
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
    mesh, output_file;
    initialize_settings = (;),
    output_settings = default_output_settings(objective),
    verbose::Bool = true
)
    pp = initialize_postprocessor(mesh, output_file, objective, output_settings)
    initialize!(objective, u, p; initialize_settings...)
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

function postprocess!(pp::PostProcessor, output_settings, n, ::AbstractSolidMechanicsObjective, u, p)
    write_times(pp, n, FEC.current_time(p.times))
    dim = size(p.field, 1)
    if output_settings.acceleration
        write_field(pp, n, _vector_field_names("acceleration", dim), p.field)
    end
    if output_settings.displacement
        write_field(pp, n, _vector_field_names("displ", dim), p.field)
    end
    if output_settings.velocity
        write_field(pp, n, _vector_field_names("velocity", dim), p.field)
    end
end

function setup_assembler(::Type{<:AbstractSolidMechanicsObjective}, mesh)
    fspace = FunctionSpace(mesh, H1Field, Lagrange)
    u = VectorFunction(fspace, "displ")
    dof = DofManager(u)
    assembler = SparseMatrixAssembler(dof; use_inplace_methods = true)
    return assembler
end

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

# helpers
function _vector_field_names(base_name::String, dim::Int)
    if dim == 2
        field_names = map(x -> "$(base_name)_$x", ("x", "y"))
    elseif dim == 3
        field_names = map(x -> "$(base_name)_$x", ("x", "y", "z"))
    else
        @assert false
    end
    return field_names
end

# logging hooks
function _step_begin_banner(::AbstractObjective, p; verbose::Bool = true)
    if verbose
        time_curr = FiniteElementContainers.current_time(p.times)
        time_start = sum(p.times.time_start)
        time_end = sum(p.times.time_end)
        str = "\n" * repeat('=', 132) * "\n"
        str = str * "Start time       = $time_start\n"
        str = str * "Current time     = $time_curr\n"
        str = str * "End time         = $time_end\n"
        str = str * "Percent complete = $(time_curr / time_end * 100)%\n"
        str = str * repeat('=', 132) * "\n"
        @info str
    end
    return nothing
end

function _step_end_banner(::AbstractObjective, p; verbose::Bool = true)
end


# include("ContactObjective.jl")
# include("ConstrainedObjective.jl")
# include("DesignObjective.jl")
include("EigenObjective.jl")
include("ExplicitDynamicsObjective.jl")
include("HeatConductionObjective.jl")
include("ImplicitDynamicsObjective.jl")
include("QuasiStaticObjective.jl")

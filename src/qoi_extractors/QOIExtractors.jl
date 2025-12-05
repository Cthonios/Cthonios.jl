# TODO
# plans for QOI extractors
# 1. have a base cache type to handle whether
# or not we need AD which is the main
# source of allocations for these structs
# 2. have input helper methods for what type of 
# calc is needed, how to reduce, and what types
# to pre-allocate for storage

abstract type AbstractQOIExtractor end

# function value end
# function value! end

# include("ScalarQOIExtractor.jl")

struct QOIExtractor{
    O          <: AbstractObjectiveCache,
    EvalFunc   <: Function,
    Reduction1 <: Function,
    Reduction2 <: Function,
    Storage,
    dSolution,
    dParameters
} <: AbstractQOIExtractor
    objective_cache::O
    func::EvalFunc
    reduction_1::Reduction1
    reduction_2::Reduction2
    storage::Storage
    dstorage::Storage # TODO allow for no gradient as well
    dU::dSolution
    dp::dParameters
end

function QOIExtractor(
    objective_cache, 
    func,
    reduction_1,
    cache_arrtype,
    cache_eltype;
    component_extractor = nothing,
    reduction_2 = identity
)
    if cache_arrtype isa Type{<:H1Field}
        @info "Requesting QOI output on reduction of H1Field"
        @assert cache_eltype isa Type{<:Number}
        storage = FiniteElementContainers.create_field(assembler(objective_cache))
    elseif cache_arrtype isa Type{<:L2QuadratureField}
        @info "Requesting QOI output on reduction of L2QuadratureField"
        if cache_eltype isa Type{<:Number}
            storage = Matrix{cache_eltype}[]
        else
            storage = StructArray{cache_eltype}[]
        end

        fspace = FiniteElementContainers.function_space(assembler(objective_cache).dof)
        for (conn, ref_fe) in zip(fspace.elem_conns, fspace.ref_fes)
            nq = num_quadrature_points(ref_fe)
            ne = size(conn, 2)
            if cache_eltype isa Type{<:Number}
                temp_output = Matrix{cache_eltype}(undef, nq, ne)
            else
                temp_output = StructArray{cache_eltype}(undef, nq, ne)
            end
            push!(storage, temp_output)
        end
    else
        @assert false "Unsupported cache type"
    end

    # TODO it won't always be a mat func
    if component_extractor !== nothing
        function mat_func(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
            val = general_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
            return val[component_extractor...]
        end
    else
        mat_func = (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) -> 
            general_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
    end

    return QOIExtractor(
        objective_cache, 
        mat_func, reduction_1, reduction_2,
        storage, deepcopy(storage),
        deepcopy(objective_cache.solution),
        deepcopy(objective_cache.parameters)
    )
end

# TODO eventually cache derivatives in qoi
function _gradient_and_value_init!(::Enzyme.ReverseMode, qoi::QOIExtractor)
    for val in values(qoi.dstorage)
        fill!(val, zero(eltype(val)))
    end
    dU = make_zero(qoi.objective_cache.solution)
    dp = make_zero(qoi.objective_cache.parameters)
    fill!(dU, zero(eltype(dU)))
    Enzyme.make_zero!(dp)
    return dU, dp
end

function _gradient_and_value!(f, df, qoi::QOIExtractor)
    Enzyme.autodiff(
        Reverse,
        _value!,
        Duplicated(f, df),
        Duplicated(qoi.storage, qoi.dstorage),
        Const(qoi.objective_cache.assembler),
        Const(qoi.func),
        Duplicated(qoi.objective_cache.solution, qoi.dU),
        Duplicated(qoi.objective_cache.parameters, qoi.dp),
        Const(qoi.reduction_1),
        Const(qoi.reduction_2)
    )
end

function gradient_props_and_value(qoi::QOIExtractor)
    f = zeros(1)
    df = ones(1) # seeding for reverse AD
    _gradient_and_value_init!(Reverse, qoi)
    _gradient_and_value!(f, df, qoi)
    return f, qoi.dp.properties
end

function gradient_u_and_value(qoi::QOIExtractor)
    f = zeros(1)
    df = ones(1) # seeding for reverse AD
    _gradient_and_value_init!(Reverse, qoi)
    _gradient_and_value!(f, df, qoi)
    return f, qoi.dU
end

function gradient_x_and_value(qoi::QOIExtractor)
    f = zeros(1)
    df = ones(1) # seeding for reverse AD
    @time _gradient_and_value_init!(Reverse, qoi)
    @time _gradient_and_value!(f, df, qoi)
    return f, qoi.dp.h1_coords
end

# TODO will need to specialize for GPU
function _value!(
    f, storage, asm, func, U, p, reduction_1, reduction_2
)
    FiniteElementContainers.assemble_quadrature_quantity!(
        storage, asm.dof, func, U, p
    )
    f_temp = mapreduce(x -> reduce(reduction_1, x), reduction_2, values(storage))
    fill!(f, f_temp)
    return nothing
end

function value(qoi::QOIExtractor)
    f = zeros(1)
    asm = assembler(qoi.objective_cache)
    U = qoi.objective_cache.solution
    p = qoi.objective_cache.parameters
    _value!(f, qoi.storage, asm, qoi.func, U, p, qoi.reduction_1, qoi.reduction_2)
    return sum(f)
end

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
    Storage
} <: AbstractQOIExtractor
    objective_cache::O
    func::EvalFunc
    reduction_1::Reduction1
    reduction_2::Reduction2
    storage::Storage
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
        function mat_func(x1, x2, x3, x4, x5, x6, x7, x8)
            val, state = general_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8)
            return val[component_extractor...], state
        end
    else
        mat_func = (x1, x2, x3, x4, x5, x6, x7, x8) -> 
            general_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8)
    end

    return QOIExtractor(
        objective_cache, 
        mat_func, reduction_1, reduction_2,
        storage
    )
end

function value_from_quadrature_field!(
    f, storage, asm, func, U, p, reduction_1, reduction_2
)
    FiniteElementContainers.assemble_quadrature_quantity!(
        storage, asm.dof, func, U, p
    )
    fill!(f, reduction_2(map(reduction_1, storage)))
end

function value(qoi::QOIExtractor)
    f = zeros(1)
    asm = assembler(qoi.objective_cache)
    U = qoi.objective_cache.solution
    p = qoi.objective_cache.parameters
    value_from_quadrature_field!(
        f, qoi.storage, asm, 
        qoi.func, U, p, 
        qoi.reduction_1, qoi.reduction_2
    )
    return sum(f)
end
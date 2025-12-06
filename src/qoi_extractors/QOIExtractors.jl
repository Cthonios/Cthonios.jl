# TODO
# plans for QOI extractors
# 1. have a base cache type to handle whether
# or not we need AD which is the main
# source of allocations for these structs
# 2. have input helper methods for what type of 
# calc is needed, how to reduce, and what types
# to pre-allocate for storage

abstract type AbstractQOIExtractor end

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
        function mat_func(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
            val = general_integrated_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
            return val[component_extractor...]
        end
    else
        mat_func = (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) -> 
            general_integrated_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
    end

    return QOIExtractor(
        objective_cache, 
        mat_func, reduction_1, reduction_2,
        storage
    )
end

# TODO will need to specialize for GPU
function _value!(
    f, storage, asm, func, U, p, reduction_1, reduction_2
)
    # TODO
    # ideally below line is not necessary
    # but we need to handle the differences between
    # condensed or not in update_for_assembly!
    FiniteElementContainers.update_field_dirichlet_bcs!(U, p.dirichlet_bcs)
    FiniteElementContainers.assemble_quadrature_quantity!(
        storage, asm.dof, func, U, p
    )
    f_temp = mapreduce(x -> reduce(reduction_1, x), reduction_2, values(storage))
    fill!(f, f_temp)
    return nothing
end

function value(qoi::QOIExtractor, U, p)
    f = zeros(1)
    asm = assembler(qoi.objective_cache)
    _value!(f, qoi.storage, asm, qoi.func, U, p, qoi.reduction_1, qoi.reduction_2)
    return f[1]
end

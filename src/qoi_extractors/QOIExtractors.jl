abstract type AbstractQOIExtractor end

struct QOIExtractor{
    FieldQOI, 
    MatQOI,
    A          <: FiniteElementContainers.AbstractAssembler,
    EvalFunc   <: Function,
    Reduction1 <: Function,
    Reduction2 <: Function,
    Storage
} <: AbstractQOIExtractor
    assembler::A
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
    is_field_qoi = false,
    is_material_qoi = false,
    reduction_2 = identity
)
    @assert is_field_qoi || is_material_qoi

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
    new_func = nothing

    if is_field_qoi
        if component_extractor !== nothing
            @assert false "Finish me"
        else
            # field_func = (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
            new_func = func
        end
    end

    if is_material_qoi
        if component_extractor !== nothing
            function mat_func(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
                val = general_integrated_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
                return val[component_extractor...]
            end
        else
            mat_func = (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) -> 
                general_integrated_material_qoi(func, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
        end

        new_func = mat_func
    end
    
    return QOIExtractor{
        is_field_qoi, is_material_qoi,
        typeof(objective_cache.assembler),
        typeof(new_func), typeof(reduction_1), typeof(reduction_2),
        typeof(storage)
    }(
        objective_cache.assembler, 
        new_func, reduction_1, reduction_2,
        storage
    )
end

# TODO need to add component extractors for _perform_reductions! methods
function _perform_reductions!(
    f, 
    qoi::QOIExtractor{false, true, A, F, R1, R2, S},
    storage
) where {A, F, R1, R2, S}
    f_temp = mapreduce(x -> reduce(qoi.reduction_1, x), qoi.reduction_2, values(storage))
    fill!(f, f_temp)
    return nothing
end

function _perform_reductions!(
    f, 
    qoi::QOIExtractor{true, false, A, F, R1, R2, S},
    storage
) where {A, F, R1, R2, S}
    # f_temp = mapreduce(x -> reduce(qoi.reduction_1, x), qoi.reduction_2, values(qoi.storage))
    # fill!(f, f_temp)
    # TODO do something useful here
    # currently not actually do any reductions
    copyto!(f.data, storage.data)
    return nothing
end

function _update_storage!(
    asm_method,
    storage, dof, func, U, p
)
    asm_method(storage, dof, func, U, p)
    return nothing
end

# deprecate
function _update_storage!(
    qoi::QOIExtractor{false, true, A, F, R1, R2, S}, 
    U, p
) where {A, F, R1, R2, S}
    _update_storage!(
        FiniteElementContainers.assemble_quadrature_quantity!,
        qoi.storage, qoi.assembler.dof, qoi.func, U, p
    )
    return nothing
end

# deprecate
function _update_storage!(
    qoi::QOIExtractor{true, false, A, F, R1, R2, S}, 
    U, p
) where {A, F, R1, R2, S}
    _update_storage!(
        FiniteElementContainers.assemble_vector!,
        qoi.storage, qoi.assembler.dof, qoi.func, U, p
    )
    return nothing
end

function _update_storage!(
    storage,
    qoi::QOIExtractor{false, true, A, F, R1, R2, S}, 
    U, p
) where {A, F, R1, R2, S}
    _update_storage!(
        FiniteElementContainers.assemble_quadrature_quantity!,
        storage, qoi.assembler.dof, qoi.func, U, p
    )
    return nothing
end

function _update_storage!(
    storage,
    qoi::QOIExtractor{true, false, A, F, R1, R2, S}, 
    U, p
) where {A, F, R1, R2, S}
    _update_storage!(
        FiniteElementContainers.assemble_vector!,
        storage, qoi.assembler.dof, qoi.func, U, p
    )
    return nothing
end

function value!(f, qoi::QOIExtractor, U, p)
    _update_storage!(qoi, U, p)
    _perform_reductions!(f, qoi, qoi.storage)
    return nothing
end

# for use in design objective
function value!(f, storage, qoi::QOIExtractor, U, p)
    _update_storage!(storage, qoi, U, p)
    _perform_reductions!(f, qoi, storage)
    return nothing
end

# TODO this out of place method will currently only
# work for scalar reductions
function value(qoi::QOIExtractor, U, p)
    f = zeros(1)
    value!(f, qoi, U, p)
    return f[1]
end

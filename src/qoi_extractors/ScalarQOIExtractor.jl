struct ScalarQOIExtractor{
    O              <: AbstractObjectiveCache,
    EvalFunc       <: Function,
    ReductionFunc1 <: Function,
    ReductionFunc2 <: Function,
    Storage,
    U,
    P
    # MS
} <: AbstractQOIExtractor
    objective_cache::O
    func::EvalFunc
    reduction_1::ReductionFunc1
    reduction_2::ReductionFunc2
    storage::Storage
    dstorage::Storage
    dU::U
    dp::P
    # map_scratch::MS
    # dmap_sratch::MS
end

function ScalarQOIExtractor(
    objective_cache, func, reduction_1,
    reduction_2=identity
)
    vals = objective_cache.assembler.scalar_quadrature_storage
    storage_vals = map(similar, values(vals))
    storage = NamedTuple{keys(vals)}(storage_vals)
    dstorage = make_zero(storage)
    dU = copy(objective_cache.solution)
    dp = deepcopy(objective_cache.parameters)
    return ScalarQOIExtractor(
        objective_cache, func, reduction_1, reduction_2,
        storage, dstorage, dU, dp
        # zeros(length(storage)), zeros(length(storage))
    )
end

# function gradient_props_and_value(qoi, U, p)
#     f = zeros(1)
#     # dfdprops = deepcopy(p.properties)
#     dfdprops = gradient_props_and_value!(f, qoi, U, p)
#     return f[1], dfdprops
# end

# function gradient_props_and_value!(
#     f, qoi, U, p
# )
#     df = make_zero(f)
#     fill!(df, one(eltype(df)))

#     dqoi = make_zero(qoi)
#     dU = make_zero(U)
#     dp = make_zero(p)
    
#     _gradient_and_value!(f, df, qoi, dqoi, U, dU, p, dp)
#     # copyto!(dfdprops, dp.properties)
#     # return dp.properties
#     return dp.times
# end

function gradient_props_and_value(qoi)
    f = zeros(1)
    df = ones(1)
    _gradient_and_value_init!(qoi)
    _gradient_and_value_eval!(f, df, qoi)
    f, qoi.dp.properties
end

function gradient_x_and_value(qoi)
    f = zeros(1)
    df = ones(1)
    _gradient_and_value_init!(qoi)
    _gradient_and_value_eval!(f, df, qoi)
    f, qoi.dp.h1_coords
end

# helper method that cover most derivatives
function _gradient_and_value_eval!(f, df, qoi)
    @time autodiff(
        Reverse,
        value!,
        Duplicated(f, df),
        Duplicated(qoi.storage, qoi.dstorage),
        Const(qoi.objective_cache.assembler),
        Const(qoi.func),
        Duplicated(qoi.objective_cache.solution, qoi.dU),
        Duplicated(qoi.objective_cache.parameters, qoi.dp),
        Const(qoi.reduction_1),
        Const(qoi.reduction_2)
        # Duplicated(qoi.map_scratch, qoi.dmap_sratch)
    )
    return nothing
end

function _gradient_and_value_init!(qoi::ScalarQOIExtractor)
    # make_zero!(qoi.dstorage)
    for val in values(qoi.dstorage)
        fill!(val, zero(eltype(val)))
    end
    # fill!(qoi.dstorage.data, zero(eltype(qoi.dstorage.data)))
    fill!(qoi.dU, zero(eltype(qoi.dU)))
    @time make_zero!(qoi.dp)
    return nothing
end

function value!(
    f, storage, asm, func, U, p, reduction_1, reduction_2
    # map_scratch
)
    FiniteElementContainers.assemble_quadrature_quantity!(
        storage, asm.dof, func, U, p
    )
    fill!(f, reduction_2(map(reduction_1, storage)))
    # map!(reduction_1, map_scratch, storage)
    # fill!(f, reduction_2(map_scratch))
    # fill!(f, reduction_2)
    return nothing
end

function value(qoi::ScalarQOIExtractor)
    f = zeros(1)
    asm = qoi.objective_cache.assembler
    U = qoi.objective_cache.solution
    p = qoi.objective_cache.parameters
    value!(
        f, qoi.storage, asm, 
        qoi.func, U, p, 
        qoi.reduction_1, qoi.reduction_2
        # qoi.map_scratch
    )
    return sum(f)
end

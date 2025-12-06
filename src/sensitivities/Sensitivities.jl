struct Sensitivity{
    Q          <: QOIExtractor,
    dStorage,
    dSolution,
    dParameters,
    Parameters <: AbstractArray,
    Solutions  <: AbstractArray{<:AbstractArray, 1}
}
    # qoi we're taking sensitivity to
    qoi::Q
    # scratch arrays for enzyme
    dstorage::dStorage
    dU::dSolution
    dp::dParameters
    # stored states
    stored_parameters::Parameters
    stored_solutions::Solutions
    # TODO need to also store state variables
    # should probably put a switch for
    # whether or not we have path-dependence
end

function Sensitivity(qoi::QOIExtractor, U, p)
    dstorage = deepcopy(qoi.storage)
    dU = deepcopy(U)
    dp = deepcopy(p)

    stored_parameters = typeof(p)[]
    stored_solutions = typeof(U)[]
    return Sensitivity(
        qoi, 
        dstorage, dU, dp,
        stored_parameters, stored_solutions
    )
end

function _gradient_and_value_init!(sens::Sensitivity, ::Enzyme.ReverseMode)
    for val in values(sens.dstorage)
        fill!(val, zero(eltype(val)))
    end
    fill!(sens.dU, zero(eltype(sens.dU)))
    Enzyme.make_zero!(sens.dp)
    return nothing
end

function _gradient_and_value!(f, df, sens::Sensitivity, U, p, ::Enzyme.ReverseMode)
    Enzyme.autodiff(
        Reverse,
        _value!,
        Duplicated(f, df),
        Duplicated(sens.qoi.storage, sens.dstorage),
        Const(sens.qoi.objective_cache.assembler),
        Const(sens.qoi.func),
        Duplicated(U, sens.dU),
        Duplicated(p, sens.dp),
        Const(sens.qoi.reduction_1),
        Const(sens.qoi.reduction_2)
    )
end

function gradient_props_and_value(sens::Sensitivity, U, p)
    f = zeros(1)
    df = ones(1)
    _gradient_and_value_init!(sens, Reverse)
    _gradient_and_value!(f, df, sens, U, p, Reverse)
    return f[1], sens.dp.properties
end

function gradient_x_and_value(sens::Sensitivity, U, p)
    f = zeros(1)
    df = ones(1)
    _gradient_and_value_init!(sens, Reverse)
    _gradient_and_value!(f, df, sens, U, p, Reverse)
    return f[1], sens.dp.h1_coords
end

function gradient_u_and_value(sens::Sensitivity, U, p)
    f = zeros(1)
    df = ones(1)
    _gradient_and_value_init!(sens, Reverse)
    _gradient_and_value!(f, df, sens, U, p, Reverse)
    return f[1], sens.dU
end

function value(sens::Sensitivity, U, p)
    return value(sens.qoi, U, p)
end

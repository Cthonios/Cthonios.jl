struct DesignObjective{
    F <: Function,
    Q,
    P,
    S
}
    func::F
    qois::Q
    dqois::Q
    stored_parameters::P
    dstored_parameters::P
    stored_solutions::S
    dstored_solutions::S
end
function DesignObjective(func, qois, U, p)
    dqois = Enzyme.make_zero(qois)
    stored_parameters = typeof(p)[]
    dstored_parameters = typeof(p)[]
    stored_solutions = typeof(U)[]
    dstored_solutions = typeof(U)[]
    return DesignObjective(
        func, qois, dqois,
        stored_parameters, dstored_parameters,
        stored_solutions, dstored_solutions
    )
end

function _gradient_and_value_init!(obj::DesignObjective, ::Enzyme.ReverseMode)
    Enzyme.make_zero!(obj.dqois)

    for val in obj.stored_parameters
        push!(obj.dstored_parameters, Enzyme.make_zero(val))
    end

    for val in obj.stored_solutions
        push!(obj.dstored_solutions, Enzyme.make_zero(val))
    end
    return nothing
end

function _gradient_and_value!(
    design_obj,
    f, df,
    qois, dqois,
    Us, dUs,
    ps, dps,
    params, dparams
)
    autodiff(
        Reverse,
        design_obj,
        Duplicated(f, df),
        Duplicated(qois, dqois),
        Duplicated(Us, dUs),
        Duplicated(ps, dps),
        Duplicated(params, dparams)
    )
end

function gradient_and_value(obj::DesignObjective, params)
    f = zeros(1)
    df = ones(1)
    dparams = Enzyme.make_zero(params)
    _gradient_and_value_init!(obj, Reverse)
    _gradient_and_value!(
        obj.func, f, df, 
        obj.qois, obj.dqois,
        obj.stored_solutions, obj.dstored_solutions,
        obj.stored_parameters, obj.dstored_parameters,
        params, dparams
    )
    return f[1], dparams
end

function value!(f, obj::DesignObjective, params)
    obj.func(f, obj.qois, obj.stored_solutions, obj.stored_parameters, params)
    return nothing
end

function value(obj::DesignObjective, params)
    f = zeros(1)
    value!(f, obj, params)
    return f[1]
end

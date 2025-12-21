struct DesignObjective{
    F <: Function,
    Q <: AbstractArray,
    P,
    S <: AbstractArray{<:AbstractField, 1},
    QS
}
    func::F
    qois::Q
    dqois::Q
    stored_parameters::P
    dstored_parameters::P
    stored_solutions::S
    dstored_solutions::S
    qoi_storages::QS
    dqoi_storages::QS
end

function DesignObjective(func, qois, U, p)
    dqois = Enzyme.make_zero(qois)
    stored_parameters = typeof(p)[]
    dstored_parameters = typeof(p)[]
    stored_solutions = typeof(U)[]
    dstored_solutions = typeof(U)[]

    # if this works, we'll likely need to reconcile
    # how to deal with multiple qois
    qoi_storages = typeof(qois[1].storage)[]
    dqoi_storages = typeof(qois[1].storage)[]

    @show typeof(qoi_storages)
    return DesignObjective(
        func, qois, dqois,
        stored_parameters, dstored_parameters,
        stored_solutions, dstored_solutions,
        qoi_storages, dqoi_storages
    )
end

function forward_problem!(obj::DesignObjective, solver, objective_cache, U, p)
    empty!(obj.stored_parameters)
    empty!(obj.dstored_parameters)
    empty!(obj.stored_solutions)
    empty!(obj.dstored_solutions)
    empty!(obj.qoi_storages)
    empty!(obj.dqoi_storages)
    initialize!(objective_cache, U, p)

    time_end = sum(p.times.time_end)

    while FiniteElementContainers.current_time(p.times) < time_end - 1e3 * eps(time_end)
        step!(solver, objective_cache, U, p; verbose=true)
        # might need to reset BCs here
        FiniteElementContainers.update_field_dirichlet_bcs!(U, p.dirichlet_bcs)

        temp = similar(U)
        # copyto!(temp.data, U.data)
        temp.data .= U.data
        # push!(obj.stored_solutions, deepcopy(U))
        push!(obj.stored_solutions, temp)
        push!(obj.stored_parameters, deepcopy(p))
    end

    # now resize stored qoi storages
    # for val in obj.qoi_storages
    # TODO currently will only work with one QOI
    for _ in 1:length(obj.stored_solutions)
        push!(obj.qoi_storages, Enzyme.make_zero(obj.qois[1].storage))
        push!(obj.dqoi_storages, Enzyme.make_zero(obj.qois[1].storage))
    end
    return nothing
end

function _gradient_and_value_init!(obj::DesignObjective, ::Enzyme.ReverseMode)
    Enzyme.make_zero!(obj.dqois)

    for val in obj.stored_parameters
        push!(obj.dstored_parameters, Enzyme.make_zero(val))
    end

    for val in obj.stored_solutions
        push!(obj.dstored_solutions, Enzyme.make_zero(val))
    end
    
    for val in obj.dqoi_storages
        Enzyme.make_zero!(val)
    end
    return nothing
end

function _gradient_and_value!(
    design_obj,
    f, df,
    qoi_storages, dqoi_storages,
    qois,
    Us, dUs,
    ps, dps,
    params, dparams
)
    autodiff(
        Reverse,
        design_obj,
        Duplicated(f, df),
        Duplicated(qoi_storages, dqoi_storages),
        Const(qois),
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
        obj.qoi_storages, obj.dqoi_storages,
        obj.qois,
        obj.stored_solutions, obj.dstored_solutions,
        obj.stored_parameters, obj.dstored_parameters,
        params, dparams
    )
    return f[1], dparams
end

function value!(f, obj::DesignObjective, params)
    obj.func(f, obj.qoi_storages, obj.qois, obj.stored_solutions, obj.stored_parameters, params)
    return nothing
end

function value(obj::DesignObjective, params)
    f = zeros(1)
    value!(f, obj, params)
    return f[1]
end

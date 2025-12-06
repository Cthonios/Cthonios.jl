struct DesignObjective{
    F <: Function,
    # Q <: QOIExtractor
    S
} <: AbstractObjective{F}
    func::F
    sens::S
end

# expected signature to be
# function (obj::DesignObjective, p)
# where p is the set of parameters we're optimizing
#
# perhapse later it shoulde
# function (obj::DesignObjective, ps)
# or
# function (obj::DesignObjective, p, args...)

# func(params, state)


function value(obj::DesignObjective, params)
    cache = obj.sens.qoi.objective_cache
    Us = obj.sens.stored_solutions
    ps = obj.sens.stored_parameters
    return obj.func(params, Us, ps)
end
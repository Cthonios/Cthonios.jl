struct ConstrainedObjective{
    C,
    O <: AbstractObjective
} <: AbstractObjective{O}
    constraint_func::C
    objective_func::O
end

struct ConstrainedObjectiveCache{
    A, O, RT, RV,
    ObjCache <: AbstractObjectiveCache{A, O, RT, RV},
    F <: Function, CP
} <: AbstractObjectiveCache{A, O, RT, RV}
    objective_cache::ObjCache
    constraint_func::F
    contact_pairs::CP
end

function _penalty(c, λ, κ)
    return ifelse(
        λ >= κ * c,
        -c * λ + 0.5 * κ * c * c,
        -0.5 * λ * λ / κ
    )
end

function _constraint_value!(
    val, contact_caches, 
    constraint_func, u, p, λ, κ
)
    c = constraint_func(contact_caches, u, p)
    # c .= _penalty()
end

function value(objective_cache, u, p, λ, κ)
    c = objective_cache.constraint_func(objective_cache.contact_pairs, u, p)
    c .= _penalty.(c, λ, κ)
    return value(objective_cache.objective_cache, u, p) + sum(c)
end

# function gradient(objective_cache, u, p, λ, κ)
#     return gradient(objective_cache.objective_cache, u, p)
# end

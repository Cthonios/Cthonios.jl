abstract type AbstractConstraint{F <: Function} end
abstract type AbstractConstraintCache{
    O <: AbstractObjectiveCache
} end

struct ConstrainedObjective{
    O <: AbstractObjective,
    C <: AbstractConstraint
}
    objective::O
    constraint::C
end

struct ConstrainedObjectiveCache{
    C <: AbstractConstraint
    O <: AbstractObjectiveCache
} <: AbstractConstraintCache{O}
    constraint::C
    objective_cache::O
end

function ConstrainedObjectiveCache(sim)

end

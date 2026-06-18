function _create_properties(physics, props)
    @assert length(physics) == length(props)
    for (k1, k2) in zip(keys(physics), keys(props))
        @assert k1 == k2
    end
    props = map(
        (x, y) -> FiniteElementContainers.create_properties(x, y), 
        values(physics), values(props)
    )
    props = map(Array, props)
    props = NamedTuple{keys(physics)}(props)
    return props
end

abstract type AbstractSimulation end

function get_objective(sim::AbstractSimulation)
    return sim.objective
end

function get_unknowns_and_parameters(sim::AbstractSimulation)
    return sim.u, sim.p
end

include("SingleDomainSimulation.jl")

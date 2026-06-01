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

function initialize!(sim::AbstractSimulation; kwargs...)
    initialize!(sim.objective, sim.u, sim.p; kwargs...)
    return nothing
end

function run!(
    sim::AbstractSimulation, solver;
    output_exodus_every = 1
)
    pp = PostProcessor(sim)
    post_process(pp, sim, 1)

    initialize!(sim)

    time_end = sim.p.times.time_end

    n = 2
    try
        while FiniteElementContainers.current_time(sim.p.times) < time_end - 1e3 * eps(time_end)
            step!(solver, sim)
            if n % output_exodus_every == 0
                post_process(pp, sim, n)
            end
            n = n + 1
        end
    finally
        close(pp)
    end
    return nothing
end

function step!(solver, sim::AbstractSimulation)
    step!(solver, sim.objective, sim.u, sim.p; verbose = solver.settings.verbose)
    return nothing
end

include("SingleDomainSimulation.jl")

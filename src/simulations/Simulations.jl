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

function initialize!(sim::AbstractSimulation; kwargs...)
    initialize!(sim.objective, sim.u, sim.p; kwargs...)
    return nothing
end

function post_process_setup(sim::AbstractSimulation, output_settings::OutputSettings)
    mesh = UnstructuredMesh(sim.mesh_file)
    V = FiniteElementContainers.function_space(assembler(sim.objective))
    nodal_vars = FEC.AbstractFunction[]
    if output_settings.acceleration
        push!(nodal_vars, VectorFunction(V, "acceleration"))
    end

    if output_settings.displacement
        push!(nodal_vars, VectorFunction(V, "displ"))
    end

    if output_settings.velocity
        push!(nodal_vars, VectorFunction(V, "velocity"))
    end

    pp = PostProcessor(mesh, sim.output_file, nodal_vars...)
    return pp
end

function post_process(pp, sim, output_settings, n)
    u_names = names(assembler(sim.objective).dof.var)
    write_times(pp, n, sim.p.times.time_current)
    if output_settings.displacement
        write_field(pp, n, u_names, sim.p.field)
    end
    # if output_settings.velocity
        # write_field(pp, n, )
end

function run!(
    sim::AbstractSimulation, solver;
    output_exodus_every = 1,
    output_settings::OutputSettings = default_output_settings(sim.objective)
)
    # pp = PostProcessor(sim)
    pp = post_process_setup(sim, output_settings)
    post_process(pp, sim, output_settings, 1)

    initialize!(sim)

    time_end = sim.p.times.time_end

    n = 2
    try
        while FiniteElementContainers.current_time(sim.p.times) < time_end - 1e3 * eps(time_end)
            step!(solver, sim)
            if n % output_exodus_every == 0
                post_process(pp, sim, output_settings, n)
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

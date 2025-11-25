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

# abstract type AbstractSimulationCache end

# function run!(sim, objective_type, solver_type)
function run!(objective_cache, solver, sim)
    pp = PostProcessor(objective_cache, sim)
    post_process(pp, objective_cache, 1)

    initialize!(objective_cache)

    p = objective_cache.parameters

    time_end = sum(p.times.time_end)

    n = 2
    try
        while FiniteElementContainers.current_time(p.times) < time_end - 1e3 * eps(time_end)
            step!(objective_cache, solver; verbose=solver.settings.verbose)
            post_process(pp, objective_cache, n)
            n = n + 1
        end
    finally
        close(pp)
    end
    return nothing
end

function _setup_simulation_common(
    sim::AbstractSimulation,
    output_file;
    # return_post_processor = true,
    use_condensed = false
)
    # handle keywords
    if output_file === nothing
        output_file = splitext(sim.mesh_file)[1] * "-output.exo"
    end

    if isfile(output_file)
        rm(output_file, force=true)
    end

    mesh = UnstructuredMesh(sim.mesh_file)
    fspace = FunctionSpace(mesh, H1Field, Lagrange)
    # fspace_q = FunctionSpace(mesh, L2QuadratureField, Lagrange)

    # check consistent field names across physics
    # and setup function space elements
    funcs = map(x -> setup_function(fspace, x), values(sim.physics))
    @assert unique(funcs) |> length == 1
    func = funcs[1]

    # setup dof manager and assembler
    dof = DofManager(func; use_condensed=use_condensed)
    assembler = SparseMatrixAssembler(dof)
    
    # finally setup parameters
    parameters = create_parameters(
        mesh,
        assembler, sim.physics, sim.properties; 
        dirichlet_bcs=sim.dirichlet_bcs,
        neumann_bcs=sim.neumann_bcs,
        times=sim.times
    )

    return assembler, parameters
end

include("SingleDomainSimulation.jl")

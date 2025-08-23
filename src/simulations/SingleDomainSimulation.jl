struct SingleDomainSimulation{M, T, P1, P2, D, N} <: AbstractSimulation
    mesh_file::M
    times::T
    physics::P1
    properties::P2
    dirichlet_bcs::D
    neumann_bcs::N
end

function SingleDomainSimulation(
    mesh_file::String,
    times::TimeStepper,
    physics::NamedTuple,
    properties::NamedTuple;
    dirichlet_bcs::Vector{<:DirichletBC} = DirichletBC[],
    neumann_bcs::Vector{<:NeumannBC} = NeumannBC[]
)
    return SingleDomainSimulation(
        mesh_file, times, physics, properties,
        dirichlet_bcs, neumann_bcs
    )
end

function simulation_cache(sim::SingleDomainSimulation)
    return SingleDomainSimulationCache(sim)
end

struct SingleDomainSimulationCache{A, P1, P2, T} <: AbstractSimulationCache{A, P1, P2, T}
    assembler::A
    parameters::P1
    post_processor::P2
    time::T
end

function SingleDomainSimulationCache(
    sim::SingleDomainSimulation;
    output_file = nothing
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

    # check consistent field names across physics
    # and setup function space elements
    funcs = map(x -> setup_function(fspace, x), values(sim.physics))
    @assert unique(funcs) |> length == 1
    func = funcs[1]

    # setup dof manager, note all dofs are unknown at init
    if length(funcs) > 1
        dof = DofManager(func...)
        post_processor = PostProcessor(mesh, output_file, func...)
    else
        dof = DofManager(func)
        post_processor = PostProcessor(mesh, output_file, func)
    end

    assembler = SparseMatrixAssembler(dof, H1Field)
    parameters = create_parameters(
        assembler, sim.physics, sim.properties; 
        dirichlet_bcs=sim.dirichlet_bcs,
        neumann_bcs=sim.neumann_bcs,
        times=sim.times
    )
    return SingleDomainSimulationCache(assembler, parameters, post_processor, TimerOutput())
end

function evolve!(sim_cache::SingleDomainSimulationCache, solver, Uu, p)
    pp = sim_cache.post_processor
    # TODO need to reset time. Probably a better way to force reset all paramters
    fill!(p.times.time_current, zero(eltype(p.times.time_current)))

    time_start = sum(p.times.time_start)
    time_end = sum(p.times.time_end)

    n = 1
    while FiniteElementContainers.current_time(p.times) < time_end
        FiniteElementContainers.update_time!(p)

        time_curr = FiniteElementContainers.current_time(p.times)
        str = "\n" * repeat('=', 132) * "\n"
        str = str * "Start time       = $time_start\n"
        str = str * "Current time     = $time_curr\n"
        str = str * "End time         = $time_end\n"
        str = str * "Percent complete = $(time_curr / time_end * 100)%\n"
        str = str * repeat('=', 132) * "\n"
        @info str

        FiniteElementContainers.update_bc_values!(p)
        Cthonios.solve!(solver, Uu, p)

        write_times(pp, n, FiniteElementContainers.current_time(p.times))
        write_field(pp, n, p.h1_field)
        n = n + 1
    end
    close(sim_cache.post_processor)
end

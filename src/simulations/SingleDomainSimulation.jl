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

struct SingleDomainSimulationCache{A, P, T} <: AbstractSimulationCache{A, P, T}
    assembler::A
    parameters::P
    time::T
end

function SingleDomainSimulationCache(sim::SingleDomainSimulation)
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
    else
        dof = DofManager(func)
    end

    assembler = SparseMatrixAssembler(dof, H1Field)
    parameters = create_parameters(
        assembler, sim.physics, sim.properties; 
        dirichlet_bcs=sim.dirichlet_bcs,
        neumann_bcs=sim.neumann_bcs,
        times=sim.times
    )
    return SingleDomainSimulationCache(assembler, parameters, TimerOutput())
end

function evolve!(sim_cache::SingleDomainSimulationCache, solver, pp, Uu, p)
    # Uu, = create_unknowns(sim_cache)
    # p = parameters(sim_cache)

    # TODO need to reset time. Probably a better way to force reset all paramters
    fill!(p.times.time_current, zero(eltype(p.times.time_current)))

    n = 1
    while FiniteElementContainers.current_time(p.times) < sum(p.times.time_end)
        FiniteElementContainers.update_time!(p)
        FiniteElementContainers.update_bc_values!(p)
        Cthonios.solve!(solver, Uu, p)

        write_times(pp, n, FiniteElementContainers.current_time(p.times))
        write_field(pp, n, p.h1_field)
        n = n + 1
    end
end

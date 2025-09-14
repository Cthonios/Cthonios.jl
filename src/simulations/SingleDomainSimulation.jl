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
    assembler, parameters, post_processor = _setup_simulation_common(sim, output_file)
    return SingleDomainSimulationCache(assembler, parameters, post_processor, TimerOutput())
end

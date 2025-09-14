abstract type AbstractSimulation end

abstract type AbstractSimulationCache{
    A  <: FiniteElementContainers.AbstractAssembler, 
    P1 <: Parameters,
    P2 <: FiniteElementContainers.PostProcessor,
    T  <: TimerOutput
} end

function FiniteElementContainers.create_unknowns(sim::AbstractSimulationCache)
    return create_unknowns(sim.assembler, H1Field)
end

function parameters(sim::AbstractSimulationCache)
    return sim.parameters
end

function run!(solver, pp)
    objective_cache = solver.objective
    initialize!(objective_cache)
    p = objective_cache.parameters
    fill!(p.times.time_current, zero(eltype(p.times.time_current)))

    time_start = sum(p.times.time_start)
    time_end = sum(p.times.time_end)

    n = 1
    while FiniteElementContainers.current_time(p.times) < time_end
        step!(objective_cache, solver)

        # post-processing
        write_times(pp, n, sum(p.times.time_current))
        write_field(pp, n, ("displ_x", "displ_y"), p.h1_field)
        n = n + 1
    end
    close(pp)
    return nothing
end

function _setup_simulation_common(
    sim::AbstractSimulation,
    output_file;
    return_post_processor = true
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

    # setup dof manager and assembler
    dof = DofManager(func)
    assembler = SparseMatrixAssembler(dof)
    
    # finally setup parameters
    parameters = create_parameters(
        mesh,
        assembler, sim.physics, sim.properties; 
        dirichlet_bcs=sim.dirichlet_bcs,
        neumann_bcs=sim.neumann_bcs,
        times=sim.times
    )

    if return_post_processor
        post_processor = PostProcessor(mesh, output_file, func)
        return assembler, parameters, post_processor
    else
        return assembler, parameters
    end
end

include("SingleDomainSimulation.jl")

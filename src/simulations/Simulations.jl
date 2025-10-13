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

# function run!(solver, pp)
function run!(sim, sim_type, solver_type)
    # mesh = 
    mesh = UnstructuredMesh(sim.mesh_file)
    V_q = FunctionSpace(mesh, L2QuadratureField, Lagrange)
    objective_cache = sim_type(sim)
    disp_var = objective_cache.assembler.dof.var

    mat_output, mat_vars = create_material_output(objective_cache, V_q, StandardMaterialOutput{Float64})
    update_material_output!(mat_output, objective_cache)

    display(mat_vars)

    pp = PostProcessor(mesh, sim.output_file, disp_var, mat_vars...)
    solver = solver_type(objective_cache)

    initialize!(objective_cache)

    p = objective_cache.parameters

    time_start = sum(p.times.time_start)
    time_end = sum(p.times.time_end)

    n = 1
    try
        while FiniteElementContainers.current_time(p.times) < time_end
            step!(objective_cache, solver)
            _post_process_common!(pp, mat_output, mat_vars, objective_cache, n)
            n = n + 1
        end
    finally
        close(pp)
    end
    return solver.timer
end

function _post_process_common!(pp, mat_output, mat_vars, objective_cache, n)
    U_names = names(assembler(objective_cache).dof.var)
    # U = objective_cache.solution
    p = objective_cache.parameters
    U = p.h1_field
    update_material_output!(mat_output, objective_cache)

    write_times(pp, n, sum(p.times.time_current))
    write_field(pp, n, U_names, U)
    # write_field(pp, n, )

    for (block, val) in pairs(mat_output)
        write_field(pp, n, String(block), "algorithmic_energy", val.algorithmic_energy)
        write_field(pp, n, String(block), "cauchy_stress", val.cauchy_stress)
        write_field(pp, n, String(block), "displacement_gradient", val.displacement_gradient)
    end
end

function _setup_simulation_common(
    sim::AbstractSimulation,
    output_file;
    return_post_processor = true,
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
    fspace_q = FunctionSpace(mesh, L2QuadratureField, Lagrange)

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

    if return_post_processor
        post_processor = PostProcessor(mesh, output_file, func)
        return assembler, parameters, post_processor
    else
        return assembler, parameters
    end
end

include("SingleDomainSimulation.jl")

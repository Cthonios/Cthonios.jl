function QuasiStaticObjective()
    return QuadratureLevelObjective(energy, residual, stiffness)
end

struct QuasiStaticObjectiveCache{
    A, O, P,
    RT, RV <: AbstractArray{RT, 1}, NF
} <: AbstractObjectiveCache2{A, O, P, RT, RV}
    assembler::A
    objective::O
    parameters::P
    #
    external_energy::RV
    internal_energy::RV
    external_force::H1Field{RT, RV, NF}
    external_force_unknowns::RV
    internal_force::H1Field{RT, RV, NF}
    internal_force_unknowns::RV
    solution::RV
    solution_old::RV
    timer::TimerOutput
end

function QuasiStaticObjectiveCache(
    # objective,
    sim
)
    objective = QuasiStaticObjective()
    assembler, parameters = _setup_simulation_common(sim, nothing; return_post_processor=false)

    RT = eltype(parameters.h1_coords)
    backend = KA.get_backend(parameters.h1_coords)
    external_energy = KA.zeros(backend, RT, 1)
    internal_energy = KA.zeros(backend, RT, 1)
    external_force = create_field(assembler)
    external_force_unknowns = create_unknowns(assembler)
    internal_force = create_field(assembler)
    internal_force_unknowns = create_unknowns(assembler)
    solution = create_unknowns(assembler)
    solution_old = create_unknowns(assembler)
    timer = TimerOutput()
    return QuasiStaticObjectiveCache(
        assembler, objective, parameters,
        external_energy, internal_energy,
        external_force, external_force_unknowns, 
        internal_force, internal_force_unknowns,
        solution, solution_old,
        timer
    )
end

# objective related methods
function gradient(o::QuasiStaticObjectiveCache, Uu, p)
    @timeit o.timer "Objective - gradient!" begin
        assemble_vector!(assembler(o), o.objective.gradient_u, Uu, p)
        copyto!(o.internal_force, assembler(o).residual_storage)
        copyto!(o.internal_force_unknowns, residual(assembler(o)))
        # assemble_vector_neumann_bc!(assembler(o), Uu, p)
        # copyto!(o.external_force, assembler(o).residual_storage .- o.internal_force)
        # copyto!(o.external_force_unknowns, residual(assembler(o)) .- o.internal_force_unknowns)
        # fill!(o.external_energy, dot(o.external_force_unknowns, solution))
    end
    return residual(assembler(o))
end

function hessian(o::QuasiStaticObjectiveCache, Uu, p)
    @timeit o.timer "Objective - hessian!" begin
        assemble_stiffness!(assembler(o), o.objective.hessian_u, Uu, p)
    end
    return stiffness(assembler(o)) |> Symmetric # needed for cholesky
end

function hvp(o::QuasiStaticObjectiveCache, Uu, p, Vu)
    @timeit o.timer "Objective - gradient!" begin
        assemble_matrix_action!(assembler(o), o.objective.hessian_u, Uu, Vu, p)
    end
    return FiniteElementContainers.hvp(assembler(o))
end

function value(o::QuasiStaticObjectiveCache, Uu, p)
    @timeit o.timer "Objective - value!" begin
        assemble_scalar!(assembler(o), o.objective.value, Uu, p)
    end
    val = mapreduce(x -> sum(x), sum, assembler(o).scalar_quadarature_storage)
    return val
end

# integration related methods
function initialize!(o::QuasiStaticObjectiveCache)
    return nothing
end

function step!(o::QuasiStaticObjectiveCache, solver)
    Uu = o.solution
    p = o.parameters
    FiniteElementContainers.update_time!(p)
    _step_begin_banner(o)
    FiniteElementContainers.update_bc_values!(p)
    solve!(solver, Uu, p)
    return nothing
end

function _step_begin_banner(o::QuasiStaticObjectiveCache)
    p = o.parameters
    time_curr = FiniteElementContainers.current_time(p.times)
    time_start = sum(p.times.time_start)
    time_end = sum(p.times.time_end)
    str = "\n" * repeat('=', 132) * "\n"
    str = str * "Start time       = $time_start\n"
    str = str * "Current time     = $time_curr\n"
    str = str * "End time         = $time_end\n"
    str = str * "Percent complete = $(time_curr / time_end * 100)%\n"
    str = str * repeat('=', 132) * "\n"
    @info str
    return nothing
end

function ImplicitDynamicsObjective()
    return QuadratureLevelObjective(energy, residual, stiffness)
end

struct ImplicitDynamicsObjectiveCache{
    A, O, P,
    RT, RV <: AbstractArray{RT, 1}, NF
} <: AbstractObjectiveCache2{A, O, P, RT, RV}
    assembler::A
    objective::O
    parameters::P
    #
    β::RT
    γ::RT
    #
    external_energy::RV
    internal_energy::RV
    kinetic_energy::RV
    external_force::H1Field{RT, RV, NF}
    external_force_unknowns::RV
    inertial_force::H1Field{RT, RV, NF}
    inertial_force_unknowns::RV
    internal_force::H1Field{RT, RV, NF}
    internal_force_unknowns::RV
    predictor_solution::RV
    predictor_solution_rate::RV
    solution::RV
    solution_rate::RV
    solution_rate_rate::RV
    solution_old::RV
    solution_rate_old::RV
    timer::TimerOutput
end

function ImplicitDynamicsObjectiveCache(
    # objective,
    sim, β=0.25, γ=0.5
)
    objective = ImplicitDynamicsObjective()
    assembler, parameters = _setup_simulation_common(sim, nothing; return_post_processor=false)

    RT = eltype(parameters.h1_coords)
    backend = KA.get_backend(parameters.h1_coords)
    external_energy = KA.zeros(backend, RT, 1)
    internal_energy = KA.zeros(backend, RT, 1)
    kinetic_energy = KA.zeros(backend, RT, 1)
    external_force = create_field(assembler)
    external_force_unknowns = create_unknowns(assembler)
    inertial_force = create_field(assembler)
    inertial_force_unknowns = create_unknowns(assembler)
    internal_force = create_field(assembler)
    internal_force_unknowns = create_unknowns(assembler)
    predictor_solution = create_unknowns(assembler)
    predictor_solution_rate = create_unknowns(assembler)
    solution = create_unknowns(assembler)
    solution_rate = create_unknowns(assembler)
    solution_rate_rate = create_unknowns(assembler)
    solution_old = create_unknowns(assembler)
    solution_rate_old = create_unknowns(assembler)
    timer = TimerOutput()
    return ImplicitDynamicsObjectiveCache(
        assembler, objective, parameters,
        β, γ,
        external_energy, internal_energy, kinetic_energy,
        external_force, external_force_unknowns, 
        inertial_force, inertial_force_unknowns,
        internal_force, internal_force_unknowns,
        predictor_solution, predictor_solution_rate,
        solution, solution_rate, solution_rate_rate,
        solution_old, solution_rate_old,
        timer
    )
end

# objective related methods
function gradient(o::ImplicitDynamicsObjectiveCache, Uu, p)
    @timeit o.timer "Objective - gradient!" begin
        assemble_vector!(assembler(o), o.objective.gradient_u, Uu, p)
        copyto!(o.internal_force, assembler(o).residual_storage)
        copyto!(o.internal_force_unknowns, residual(assembler(o)))
        # assemble_vector_neumann_bc!(assembler(o), Uu, p)
        # copyto!(o.external_force, assembler(o).residual_storage .- o.internal_force)
        # copyto!(o.external_force_unknowns, residual(assembler(o)) .- o.internal_force_unknowns)
        # fill!(o.external_energy, dot(o.external_force_unknowns, solution))

        # inertial force
        assemble_mass!(assembler(o), mass, Uu, p)
        M = FiniteElementContainers.mass(assembler(o))
        # copyto!(o.inertial_force_unknowns, M * o.solution_rate_rate)
        mul!(o.inertial_force_unknowns, M, o.solution_rate_rate) # this can't be right
    end
    return residual(assembler(o)) .+ o.inertial_force_unknowns
end

function hessian(o::ImplicitDynamicsObjectiveCache, Uu, p)
    @timeit o.timer "Objective - hessian!" begin
        assemble_mass!(assembler(o), mass, Uu, p)
        assemble_stiffness!(assembler(o), o.objective.hessian_u, Uu, p)
    end
    # return stiffness(assembler(o)) |> Symmetric # needed for cholesky
    Δt = FiniteElementContainers.time_step(parameters(o).times)
    β = o.β
    # γ = o.γ
    H = stiffness(assembler(o)) + FiniteElementContainers.mass(assembler(o)) / Δt / Δt / β
    return H |> Symmetric
end

# TODO hvp

function value(o::ImplicitDynamicsObjectiveCache, Uu, p)
    # @timeit o.timer "Objective - value!" begin
    # end

    # strain energy
    assemble_scalar!(assembler(o), o.objective.value, Uu, p)
    val_internal = mapreduce(x -> sum(x), sum, assembler(o).scalar_quadarature_storage)
    fill!(o.internal_energy, val_internal)

    # kinetic energy
    assemble_scalar!(assembler(o), kinetic_energy, o.solution_rate, p)
    val_inertial = mapreduce(x -> sum(x), sum, assembler(o).scalar_quadarature_storage)
    fill!(o.kinetic_energy, val_inertial)

    # TODO add in external energy
    return val_internal + val_inertial
end

# integration methods
function initialize!(o::ImplicitDynamicsObjectiveCache)
    # set initial energy
    value(o, o.solution, o.parameters)

    # do assembly
    # assemble_scalar!(assembler(o), kinetic_energy, o.solution_rate, o.parameters)
    assemble_mass!(assembler(o), mass, o.solution_rate, o.parameters)
    # assemble_scalar!(assembler(o), kinetic_energy, o.solution_rate, o.parameters)
    M = FiniteElementContainers.mass(assembler(o))
    f = gradient(o, o.solution, o.parameters)
    f_inertial = -f
    a = M \ f_inertial

    # set stuff
    copyto!(o.solution_rate_rate, a)
end

function step!(o::ImplicitDynamicsObjectiveCache, solver)
    Uu, p = o.solution, o.parameters
    FiniteElementContainers.update_time!(p)
    _step_begin_banner(o)
    FiniteElementContainers.update_bc_values!(p)
    _predict!(o)
    solve!(solver, Uu, p)
    _correct!(o)
    return nothing
end

function _correct!(o::ImplicitDynamicsObjectiveCache)
    # TODO
    u = o.solution
    v = o.solution_rate
    a = o.solution_rate_rate
    u_pre = o.predictor_solution
    v_pre = o.predictor_solution_rate
    Δt = FiniteElementContainers.time_step(parameters(o).times)

    copyto!(a, (u - u_pre) / o.β / Δt / Δt)
    copyto!(v, v_pre + o.γ * Δt * a)
    return nothing
end

# this updates displacement and velocity only
function _predict!(o::ImplicitDynamicsObjectiveCache)
    copyto!(o.predictor_solution, o.solution)
    copyto!(o.predictor_solution_rate, o.solution_rate)

    u = o.predictor_solution
    v = o.predictor_solution_rate
    a = o.solution_rate_rate
    Δt = FiniteElementContainers.time_step(parameters(o).times)

    # calculate predictor and store in current solution
    u .+= Δt * v + (0.5 - o.β) * Δt * Δt * a
    v .+= (1.0 - o.γ) * Δt * a
    copyto!(o.solution, u)
    copyto!(o.solution_rate, v)
    return nothing
end

function _step_begin_banner(o::ImplicitDynamicsObjectiveCache)
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

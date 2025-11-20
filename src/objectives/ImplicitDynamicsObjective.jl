struct ImplicitDynamicsObjective{
    F1 <: Function,
    F2 <: Function,
    F3 <: Function
} <: AbstractObjective{F1}
    value::F1
    gradient_u::F2
    hessian_u::F3
end

function ImplicitDynamicsObjective()
    return ImplicitDynamicsObjective(energy, residual, stiffness)
end

struct ImplicitDynamicsObjectiveCache{
    A, O, P,
    RT, RV <: AbstractArray{RT, 1}, NF
} <: AbstractObjectiveCache{A, O, P, RT, RV}
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
    inertial_force::H1Field{RT, RV, NF}
    internal_force::H1Field{RT, RV, NF}
    predictor_solution::H1Field{RT, RV, NF}
    predictor_solution_rate::H1Field{RT, RV, NF}
    solution::H1Field{RT, RV, NF}
    solution_rate::H1Field{RT, RV, NF}
    solution_rate_rate::H1Field{RT, RV, NF}
    solution_old::H1Field{RT, RV, NF}
    solution_rate_old::H1Field{RT, RV, NF}
    #
    value::RV
    gradient::H1Field{RT, RV, NF}
    #
    timer::TimerOutput
end

function ImplicitDynamicsObjectiveCache(sim, β=0.25, γ=0.5)
    objective = ImplicitDynamicsObjective()
    assembler, parameters = _setup_simulation_common(
        sim, nothing; 
        return_post_processor=false,
        use_condensed=true
    )

    RT = eltype(parameters.h1_coords)
    backend = KA.get_backend(parameters.h1_coords)

    external_energy = KA.zeros(backend, RT, 1)
    internal_energy = KA.zeros(backend, RT, 1)
    kinetic_energy = KA.zeros(backend, RT, 1)

    external_force = create_field(assembler)
    inertial_force = create_field(assembler)
    internal_force = create_field(assembler)

    predictor_solution = create_field(assembler)
    predictor_solution_rate = create_field(assembler)

    solution = create_field(assembler)
    solution_rate = create_field(assembler)
    solution_rate_rate = create_field(assembler)
    
    solution_old = create_field(assembler)
    solution_rate_old = create_field(assembler)

    value = KA.zeros(backend, RT, 1)
    gradient = create_field(assembler)

    timer = TimerOutput()

    return ImplicitDynamicsObjectiveCache(
        assembler, objective, parameters,
        β, γ,
        external_energy, internal_energy, kinetic_energy,
        external_force, inertial_force, internal_force,
        predictor_solution, predictor_solution_rate,
        solution, solution_rate, solution_rate_rate,
        solution_old, solution_rate_old,
        value, gradient,
        timer
    )
end

function gradient(o::ImplicitDynamicsObjectiveCache, U, p)
    RT = eltype(U)
    assemble_vector!(assembler(o), o.objective.gradient_u, U, p)
    # copyto!(o.internal_force, residual(assembler(o)))
    copyto!(o.internal_force, assembler(o).residual_storage)
    fill!(o.external_force, zero(RT))
    fill!(assembler(o).residual_storage, zero(RT))
    # fill!(assembler(o).residual_unknowns, zero(RT))
    assemble_vector_neumann_bc!(assembler(o), U, p)
    copyto!(o.external_force, assembler(o).residual_storage)

    Δt = FiniteElementContainers.time_step(parameters(o).times)
    β = o.β

    # assume mass is already assembled so we don't have to do it here
    M = FiniteElementContainers.mass(assembler(o))
    mul!(o.inertial_force.data, M, o.solution_rate_rate.data)
    # mul!(o.inertial_force.data, M, (o.solution.data - o.predictor_solution.data) / Δt / Δt / β)

    o.gradient .= o.internal_force - o.external_force + o.inertial_force
    return o.gradient.data    
end

function hessian(o::ImplicitDynamicsObjectiveCache, U, p)
    # also assume mass is already assembled here
    U_pred = o.predictor_solution.data
    assemble_stiffness!(assembler(o), o.objective.hessian_u, U, p)

    Δt = FiniteElementContainers.time_step(parameters(o).times)
    β = o.β

    # assemble into hessian storage
    asm = assembler(o)
    asm.hessian_storage .= asm.stiffness_storage .+ asm.mass_storage / Δt / Δt / β
    H = FiniteElementContainers.hessian(asm)
    return H
end

function initialize!(o::ImplicitDynamicsObjectiveCache)
    # set initial energy
    # value(o, o.solution, o.parameters)

    # do assembly
    # assemble_scalar!(assembler(o), kinetic_energy, o.solution_rate, o.parameters)
    assemble_mass!(assembler(o), mass, o.solution_rate, o.parameters)
    # assemble_scalar!(assembler(o), kinetic_energy, o.solution_rate, o.parameters)
    M = FiniteElementContainers.mass(assembler(o))
    f = gradient(o, o.solution, o.parameters)
    f_inertial = -f
    a = M \ f_inertial

    # set stuff
    copyto!(o.solution_rate_rate.data, a)
    # @assert false
    return nothing
end

function step!(o::ImplicitDynamicsObjectiveCache, solver)
    U, p = o.solution, o.parameters
    V, A = o.solution_rate, o.solution_rate_rate
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p)
    FiniteElementContainers.update_field_dirichlet_bcs!(U, V, A, p.dirichlet_bcs)
    _step_begin_banner(o)
    _predict!(o)
    # FiniteElementContainers.update_field_dirichlet_bcs!(U, V, A, p.dirichlet_bcs)
    # FiniteElementContainers.update_field_dirichlet_bcs!(o.predictor_solution, o.predictor_solution_rate, A, p.dirichlet_bcs)
    solve!(solver, U.data, p)
    _correct!(o)

    # copyto!(o.solution_old, o.solution)
    # copyto!(o.solution_rate_old, o.solution_rate)
    return nothing
end

function _correct!(o::ImplicitDynamicsObjectiveCache)
    # TODO
    u = o.solution
    v = o.solution_rate
    # v = o.predictor_solution_rate
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
    u .+= Δt * v .+ (0.5 - o.β) * Δt * Δt * a
    v .+= (1.0 - o.γ) * Δt * a

    FiniteElementContainers.update_field_dirichlet_bcs!(u, v, o.parameters.dirichlet_bcs)

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

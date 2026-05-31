struct ImplicitDynamicsObjective{
    F1 <: Function,
    F2 <: Function,
    F3 <: Function,
    F4 <: Function
} <: AbstractObjective{F1}
    value::F1
    gradient_u::F2
    hessian_u::F3
    hvp_u::F4
end

function ImplicitDynamicsObjective(; use_inplace_methods = true)
    if use_inplace_methods
        return ImplicitDynamicsObjective(energy, residual!, stiffness!, stiffness_action!)
    else
        return ImplicitDynamicsObjective(energy, residual, stiffness, stiffness_action)
    end
end

mutable struct ImplicitDynamicsObjectiveCache{
    A, O,
    RT, RV <: AbstractArray{RT, 1}
} <: AbstractObjectiveCache{A, O, RT, RV}
    assembler::A
    objective::O
    #
    β::RT
    γ::RT
    α_hht::RT
    #
    external_energy::RT
    internal_energy::RT
    kinetic_energy::RT
    #
    v::RV
    a::RV
    u_prev::RV
    v_prev::RV
    a_prev::RV
    u_pred::RV
    du::RV
    R_eff::RV
    F_int_n::RV
    #
    hvp_scratch::RV
    #
    timer::TimerOutput
end

function ImplicitDynamicsObjectiveCache(
    assembler, objective::ImplicitDynamicsObjective,
    β = 0.25, γ = 0.5, α_hht = 0.0
)
    RT = eltype(assembler.constraint_storage)

    external_energy = zero(RT)
    internal_energy = zero(RT)
    kinetic_energy = zero(RT)

    v = create_unknowns(assembler)
    a = create_unknowns(assembler)
    u_prev = create_unknowns(assembler)
    v_prev = create_unknowns(assembler)
    a_prev = create_unknowns(assembler)
    u_pred = create_unknowns(assembler)
    du = create_unknowns(assembler)
    R_eff = create_unknowns(assembler)
    F_int_n = create_unknowns(assembler)

    hvp_scratch = create_unknowns(assembler)

    timer = TimerOutput()

    return ImplicitDynamicsObjectiveCache(
        assembler, objective,
        β, γ, α_hht,
        external_energy, internal_energy, kinetic_energy,
        v, a,
        u_prev, v_prev, a_prev,
        u_pred, du, R_eff, F_int_n,
        hvp_scratch,
        timer
    )
end

function setup_cache(assembler, objective::ImplicitDynamicsObjective, args...)
    return ImplicitDynamicsObjectiveCache(assembler, objective, args...)
end

function initialize!(
    o::ImplicitDynamicsObjectiveCache, u, p;
    displ_ics = nothing,
    vel_ics = nothing
)
    p.times.time_current = zero(p.times.time_current)
    # Δt = _compute_stable_dt(assembler(o), p, o.CFL)
    # p.times.Δt = Δt

    # apply initial conditions 
    # TODO make this better in FEC
    if displ_ics !== nothing
        u_all = create_field(o)
        FiniteElementContainers.update_ic_values!(displ_ics, p.coords)
        FiniteElementContainers.update_field_ics!(u_all, displ_ics)
        FiniteElementContainers.extract_field_unknowns!(u, assembler(o).dof, u_all)
    end

    if vel_ics !== nothing
        v_all = create_field(o)
        FiniteElementContainers.update_ic_values!(vel_ics, p.coords)
        FiniteElementContainers.update_field_ics!(v_all, vel_ics)
        FiniteElementContainers.extract_field_unknowns!(o.v, assembler(o).dof, v_all)
    end

    # setup mass matrix once
    assemble_mass!(assembler(o), mass, u, p)
    M = mass(assembler(o))

    # rhs = -gradient(o, u, p)
    rhs = -_internal_force(o, u, p)
    a0, stats = Krylov.cg(M, rhs; atol = 0.0, rtol = 1e-12, verbose = 0)
    o.a .= a0
    return nothing
end

function step!(solver, o::ImplicitDynamicsObjectiveCache, u, p; verbose = true)
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p, assembler(solver.objective_cache))
    Δt = FiniteElementContainers.time_step(p)
    _step_begin_banner(o, p; verbose = verbose)

    # predict
    c_M = one(Δt) / (o.β * Δt * Δt)
    o.u_prev .= u
    o.v_prev .= o.v
    o.a_prev .= o.a
    @. u = o.u_prev + Δt * o.v_prev + Δt * Δt * (0.5 - o.β) * o.a_prev
    @. o.v = o.v_prev + Δt * (one(o.γ) - o.γ) * o.a_prev
    o.u_pred .= u
    fill!(o.du, zero(eltype(o.du)))

    # # calculate du
    # @. o.du = u - o.u_pred

    # # get gradient
    # R_int = _internal_force(o, u, p)

    # # need mass action
    # assemble_matrix_free_action!(assembler(o), mass_action, u, o.du, p)
    # M_du = FiniteElementContainers.hvp(assembler(o), o.du)
    # @. o.R_eff = -((1 + o.α_hht) * R_int + c_M * M_du - o.α_hht * o.F_int_n)
    # assemble_matrix_action!(assembler(o), mass_action!, u, o.du, p)

    # solve
    solve!(solver, u, p)

    # correct
    @. o.a = c_M * (u - o.u_pred)
    @. o.v += Δt * o.γ * o.a

    # update some values
    value(o, u, p)

    _step_end_banner(o, p; verbose = verbose)
    return nothing
end

function gradient(o::ImplicitDynamicsObjectiveCache, u, p)
    Δt = FiniteElementContainers.time_step(p)
    c_M = one(Δt) / (o.β * Δt * Δt)

    # get gradient
    R_int = _internal_force(o, u, p)

    # calculate du
    @. o.du = u - o.u_pred

    # need mass action
    assemble_matrix_free_action!(assembler(o), mass_action, u, o.du, p)
    M_du = FiniteElementContainers.hvp(assembler(o), o.du)
    @. o.R_eff = ((1 + o.α_hht) * R_int + c_M * M_du - o.α_hht * o.F_int_n)
    return o.R_eff
end

function hessian(o::ImplicitDynamicsObjectiveCache, u, p; symmetric = false)
    # assumes assemble_mass! has already been called
    Δt = FiniteElementContainers.time_step(p)
    c_M = one(Δt) / (o.β * Δt * Δt)
    asm = assembler(o)
    assemble_stiffness!(asm, o.objective.hessian_u, u, p)
    @. asm.stiffness_storage += c_M * asm.mass_storage

    if symmetric
        H = stiffness(asm) |> Symmetric
    else
        H = stiffness(asm)
    end
    return H
end

function hvp(o::ImplicitDynamicsObjectiveCache, u, v, p)
    Δt = FiniteElementContainers.time_step(p)
    c_M = one(Δt) / (o.β * Δt * Δt)

    assemble_matrix_free_action!(assembler(o), o.objective.hvp_u, u, v, p)
    Kv = FiniteElementContainers.hvp(assembler(o), v)

    o.hvp_scratch .= Kv

    assemble_matrix_free_action!(assembler(o), mass_action, u, v, p)

    Mv = FiniteElementContainers.hvp(assembler(o), v)

    # (K + c_M M)v
    @. o.hvp_scratch += c_M * Mv

    return o.hvp_scratch
end

# function value(o::ImplicitDynamicsObjectiveCache, u, p)
#     # TODO external energy
#     o.external_energy = zero(o.external_energy)
#     assemble_scalar!(assembler(o), o.objective.value, u, p)
#     o.internal_energy = sum(assembler(o).scalar_quadrature_storage)
#     # map!((m, v) -> m * v * v, o.value_scratch, o.lumped_mass, o.v)
#     # o.kinetic_energy = reduce(+, o.value_scratch)
#     assemble_matrix_free_action!(assembler(o), mass_action, u, o.v, p)

# end

function value(o::ImplicitDynamicsObjectiveCache, u, p)
    # TODO external energy
    o.external_energy = zero(o.external_energy)
    assemble_scalar!(assembler(o), o.objective.value, u, p)
    o.internal_energy = sum(assembler(o).scalar_quadrature_storage)
    @. o.du = u - o.u_pred
    assemble_matrix_free_action!(assembler(o), mass_action, u, o.du, p)
    Mdu = FiniteElementContainers.hvp(assembler(o), o.du)

    Δt = FiniteElementContainers.time_step(p)
    c_M = one(Δt) / (o.β * Δt * Δt)

    o.kinetic_energy = 0.5 * c_M * dot(o.du, Mdu)

    return o.external_energy + o.internal_energy + o.kinetic_energy
end

function _internal_force(o::ImplicitDynamicsObjectiveCache, u, p)
    # o.inertial_force .= lumped_mass .* o.a
    assemble_vector!(assembler(o), o.objective.gradient_u, u, p)
    # o.internal_force .= assembler(o).residual_storage
    assemble_vector_neumann_bc!(assembler(o), u, p)
    assemble_vector_source!(assembler(o), u, p)
    # o.external_force .= assembler(o).residual_storage .- o.internal_force
    return residual(assembler(o))
end

function _step_begin_banner(o::ImplicitDynamicsObjectiveCache, p; verbose = true)
    if verbose
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
    end
    return nothing
end

function _step_end_banner(o::ImplicitDynamicsObjectiveCache, p; verbose::Bool = true)
    if verbose
        @info "External energy        = $(o.external_energy)"
        @info "Internal energy        = $(o.internal_energy)"
        @info "Kinetic energy         = $(o.kinetic_energy)"
        @info "Total energy           = $(o.external_energy + o.internal_energy + o.kinetic_energy)"
    end
end
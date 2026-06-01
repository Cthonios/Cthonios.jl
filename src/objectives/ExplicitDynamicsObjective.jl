# currently uses central difference method

mutable struct ExplicitDynamicsObjective{RT, RV, A} <: AbstractSolidMechanicsObjective{RT, RV, A}
    assembler::A
    #
    γ::RT
    CFL::RT
    #
    external_energy::RT
    internal_energy::RT
    kinetic_energy::RT
    v::RV
    a::RV
    R_eff::RV
    lumped_mass::RV
    #
    value_scratch::RV
    #
    step_number::Int
    #
    timer::TimerOutput
end

function ExplicitDynamicsObjective(
    assembler, 
    CFL, γ = 0.5
)
    RT = eltype(assembler.constraint_storage)

    external_energy = zero(RT)
    internal_energy = zero(RT)
    kinetic_energy = zero(RT)

    v = create_unknowns(assembler)
    a = create_unknowns(assembler)
    R_eff = create_unknowns(assembler)
    lumped_mass = create_unknowns(assembler)
    value_scratch = create_unknowns(assembler)
    timer = TimerOutput()

    return ExplicitDynamicsObjective(
        assembler,
        γ, CFL,
        external_energy, internal_energy, kinetic_energy,
        v, a,
        R_eff, lumped_mass,
        value_scratch,
        0, timer
    )
end

function initialize!(
    o::ExplicitDynamicsObjective, u, p;
    displ_ics = nothing,
    vel_ics = nothing
)
    # fill!(p.times.time_current, zero(eltype(p.times.time_current)))
    p.times.time_current = zero(p.times.time_current)
    Δt = _compute_stable_dt(assembler(o), p, o.CFL)
    p.times.Δt = Δt

    # setup lumped mass once
    FiniteElementContainers.assemble_diagonal!(assembler(o), mass, u, p)
    o.lumped_mass .= FiniteElementContainers.diagonal(assembler(o))

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

    # calculate mechanics stuff at beginning
    value(o, u, p)

    # calculate rhs
    rhs = -gradient(o, u, p)
    # TODO external force
    # external_force = create_field(o.assembler)
    # o.external_force .= external_force
    # o.inertial_force .= internal_force .- external_force
    # o.solution_rate_rate.data .= o.inertial_force.data ./ o.lumped_mass
    m = lumped_mass(o, u, p)
    o.a .= rhs ./ m
    return nothing
end

function step!(solver, o::ExplicitDynamicsObjective, u, p; verbose = true)
    @timeit o.timer "ExplicitDynamicsObjective - step!" begin
        FiniteElementContainers.update_time!(p)
        FiniteElementContainers.update_bc_values!(p, assembler(solver.objective_cache))

        # unpack some stuff for convenience
        Δt = FiniteElementContainers.time_step(p)

        # predict
        @timeit o.timer "predict" begin
            @. u   += Δt * o.v + 0.5 * Δt^2 * o.a
            @. o.v += (one(o.γ) - o.γ) * Δt * o.a
        end

        # "solve"
        @timeit o.timer "solve" begin
            @timeit o.timer "gradient" begin
                R_eff = gradient(o, u, p)
            end
            m = lumped_mass(o, u, p)
            @timeit o.timer "solve" begin
                @. o.a = -(R_eff / m)
            end
        end

        # correct
        @timeit o.timer "correct" begin
            @. o.v += o.γ * Δt * o.a
        end

        # evaluate stuff at end of step
        value(o, u, p)
    end
    _step_log(o, p)
end

function gradient(o::ExplicitDynamicsObjective, u, p)
    # o.inertial_force .= lumped_mass .* o.a
    assemble_vector!(assembler(o), residual!, u, p)
    # o.internal_force .= assembler(o).residual_storage
    assemble_vector_neumann_bc!(assembler(o), u, p)
    assemble_vector_source!(assembler(o), u, p)
    # o.external_force .= assembler(o).residual_storage .- o.internal_force
    return residual(assembler(o))
end

function lumped_mass(o::ExplicitDynamicsObjective, u, p)
    return o.lumped_mass
end

function value(o::ExplicitDynamicsObjective, u, p)
    # TODO external energy
    o.external_energy = zero(o.external_energy)
    assemble_scalar!(assembler(o), energy, u, p)
    o.internal_energy = sum(assembler(o).scalar_quadrature_storage)
    map!((m, v) -> m * v * v, o.value_scratch, o.lumped_mass, o.v)
    o.kinetic_energy = reduce(+, o.value_scratch)
    return o.external_energy + o.internal_energy + o.kinetic_energy
end

# copied from Carina.jl
function _compute_stable_dt(asm, p, CFL)
    fspace = FiniteElementContainers.function_space(asm.dof)

    # Pre-allocate per-block storage for element char lengths (nq × nelem)
    char_len_storage = []
    for (b, ref_fe) in enumerate(fspace.ref_fes)
        nquad = num_cell_quadrature_points(ref_fe)
        nelem = num_elements(fspace, b)
        push!(char_len_storage, zeros(Float64, nquad, nelem))
    end
    char_len_storage = NamedTuple{keys(fspace.ref_fes)}(char_len_storage)

    # Assemble per-element char lengths on device
    U_zeros = zeros(Float64, length(asm.dof.unknown_dofs))
    FiniteElementContainers.assemble_quadrature_quantity!(
        char_len_storage, nothing, asm.dof,
        characteristic_element_length,
        U_zeros, p
    )

    # Min-reduction over all blocks
    stable_dt = Inf
    for (b, (block_physics, block_storage, props)) in enumerate(zip(
        values(p.physics), values(char_len_storage), values(p.properties),
    ))
        ρ   = props[1]
        M   = p_wave_modulus(block_physics.constitutive_model, props)
        c_p = sqrt(M / ρ)
        h_min = minimum(block_storage)   # GPU-native reduction if on device
        block_dt = CFL * h_min / c_p
        stable_dt = min(stable_dt, block_dt)
    end
    return stable_dt
end

function _step_log(o::ExplicitDynamicsObjective, p)
    if o.step_number % 1000 == 0
        @info "Step      Current      Time        Internal    External    Kinetic"
        @info "Number    Time         Increment   Energy      Energy      Energy"
    end

    if o.step_number % 100 == 0
        str = @sprintf "%8d  %.4e   %.4e  %.4e  %.4e  %.4e" o.step_number p.times.time_current p.times.Δt[1] o.internal_energy[1] o.external_energy[1] o.kinetic_energy[1]
        @info str
    end

    o.step_number = o.step_number + 1
end

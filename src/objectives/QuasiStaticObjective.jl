mutable struct QuasiStaticObjective{RT, RV, A, NF} <: AbstractSolidMechanicsObjective{RT, RV, A}
    assembler::A
    #
    external_energy::RT
    internal_energy::RT
    external_force::H1Field{RT, RV, NF}
    internal_force::H1Field{RT, RV, NF}
    solution_old::H1Field{RT, RV, NF}
    #
    timer::TimerOutput
end

function QuasiStaticObjective(assembler)
    RT = eltype(assembler.constraint_storage)
    external_energy = zero(RT)
    internal_energy = zero(RT)
    external_force = create_field(assembler)
    internal_force = create_field(assembler)
    solution_old = create_field(assembler)

    timer = TimerOutput()

    return QuasiStaticObjective(
        assembler,
        external_energy, internal_energy,
        external_force, internal_force,
        solution_old,
        timer
    )
end

# objective hooks
function gradient(o::QuasiStaticObjective, u, p)
    assemble_vector!(assembler(o), residual!, u, p)
    o.internal_force .= assembler(o).residual_storage
    assemble_vector_neumann_bc!(assembler(o), u, p)
    o.external_force .= assembler(o).residual_storage .- o.internal_force
    return residual(assembler(o))
end

function hessian(o::QuasiStaticObjective, u, p; symmetric = true)
    assemble_stiffness!(assembler(o), stiffness!, u, p)
    if symmetric
        H = stiffness(assembler(o)) |> Symmetric
    else
        H = stiffness(assembler(o))
    end
    return H
end

function hvp(o::QuasiStaticObjective, u, v, p)
    asm = assembler(o)
    assemble_matrix_action!(asm, stiffness_action!, u, v, p)
    return FiniteElementContainers.hvp(asm, v)
end

function mass_matrix(o::QuasiStaticObjective, u, p)
    assemble_mass!(assembler(o), mass!, u, p)
    M = FiniteElementContainers.mass(assembler(o))
    return M |> Symmetric
end

# NOTE I think this assumes o.external_force is already populated
function value(o::QuasiStaticObjective, U, p)
    assemble_scalar!(assembler(o), energy, U, p)
    o.internal_energy = reduce(+, assembler(o).scalar_quadrature_storage)
    # TODO need an actual neumann energy method in FEContainers
    # fill!(o.external_energy, dot(o.external_force, U))
    o.external_energy = dot(o.external_force, p.field)
    return o.external_energy + o.internal_energy
end

# integrator hooks
function initialize!(::QuasiStaticObjective, u, p)
    p.times.time_current = zero(typeof(p.times.time_current))
    return nothing
end

function step!(solver, o::QuasiStaticObjective, u, p; verbose = true)
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p, solver.objective_cache.assembler)
    _step_begin_banner(o, p; verbose = verbose)
    solve!(solver, u, p)

    # update values at end of step
    gradient(o, u, p)
    value(o, u, p)

    # update old solution
    copyto!(o.solution_old, p.field)
    copyto!(p.state_old, p.state_new)

    _step_end_banner(o, p; verbose = verbose)
    return nothing
end

# logging hooks
function _step_begin_banner(::QuasiStaticObjective, p; verbose::Bool = true)
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

function _step_end_banner(o::QuasiStaticObjective, p; verbose::Bool = true)
    if verbose
        @info "External energy        = $(o.external_energy)"
        @info "Internal energy        = $(o.internal_energy)"
        @info "Total energy           = $(o.external_energy + o.internal_energy)"
    end
end

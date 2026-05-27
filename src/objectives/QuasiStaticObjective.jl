struct QuasiStaticObjective{
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

function QuasiStaticObjective(; use_inplace_methods::Bool = false)
    if use_inplace_methods
        return QuasiStaticObjective(energy, residual!, stiffness!, stiffness_action!)
    else
        return QuasiStaticObjective(energy, residual, stiffness, stiffness_action)
    end
end

mutable struct QuasiStaticObjectiveCache{
    A, O,
    RT, RV <: AbstractArray{RT, 1},
    NF
} <: AbstractObjectiveCache{A, O, RT, RV}
    assembler::A
    objective::O
    #
    external_energy::RT
    internal_energy::RT
    external_force::H1Field{RT, RV, NF}
    internal_force::H1Field{RT, RV, NF}
    solution_old::H1Field{RT, RV, NF}
    # solver helpers
    value::RT
    gradient::H1Field{RT, RV, NF}
    #
    timer::TimerOutput
end

function QuasiStaticObjectiveCache(
    assembler,
    objective::QuasiStaticObjective
)
    RT = eltype(assembler.constraint_storage)
    external_energy = zero(RT)
    internal_energy = zero(RT)
    external_force = create_field(assembler)
    internal_force = create_field(assembler)
    solution_old = create_field(assembler)

    value = zero(RT)
    gradient = create_field(assembler)

    timer = TimerOutput()

    return QuasiStaticObjectiveCache(
        assembler, objective,
        external_energy, internal_energy,
        external_force, internal_force,
        solution_old,
        value, gradient,
        timer
    )
end

# cache hook
function setup_cache(assembler, objective::QuasiStaticObjective)
    return QuasiStaticObjectiveCache(assembler, objective)
end

# objective hooks
function gradient(o::QuasiStaticObjectiveCache, U, p)
    assemble_vector!(assembler(o), o.objective.gradient_u, U, p)
    o.internal_force .= assembler(o).residual_storage
    assemble_vector_neumann_bc!(assembler(o), U, p)
    o.external_force .= assembler(o).residual_storage .- o.internal_force
    o.gradient .= assembler(o).residual_storage
    return residual(assembler(o))
end

function hessian(o::QuasiStaticObjectiveCache, U, p; symmetric = true)
    assemble_stiffness!(assembler(o), o.objective.hessian_u, U, p)
    if symmetric
        H = stiffness(assembler(o)) |> Symmetric
    else
        H = stiffness(assembler(o))
    end
    return H
end

function hvp(o::QuasiStaticObjectiveCache, U, V, p)
    # assemble_matrix_action!(assembler(o), o.objective.hessian_u, U, V, p)
    asm = assembler(o)
    if FiniteElementContainers._use_inplace_methods(asm)
        assemble_matrix_action!(asm, o.objective.hvp_u, U, V, p)
    else
        assemble_matrix_free_action!(asm, o.objective.hvp_u, U, V, p)
    end
    return FiniteElementContainers.hvp(asm, V)
end

# TODO review this
function mass_matrix(o::QuasiStaticObjectiveCache, U, p)
    assemble_mass!(assembler(o), o.objective.hessian_u, U, p)
    # M = FiniteElementContainers.mass(assembler(o)) |> Symmetric
    # really passing displacement below instead of velocity to save
    # on having to store the solution_rate for quasi-statics
    # this method is mainly so you can construct a mass matrix
    # if you want one. Should add abstract types in the future
    assemble_mass!(assembler(o), mass, o.solution, o.parameters)
    M = FiniteElementContainers.mass(assembler(o))
    return M |> Symmetric
end

# NOTE I think this assumes o.external_force is already populated
function value(o::QuasiStaticObjectiveCache, U, p)
    assemble_scalar!(assembler(o), o.objective.value, U, p)
    o.internal_energy = reduce(+, assembler(o).scalar_quadrature_storage)
    # TODO need an actual neumann energy method in FEContainers
    # fill!(o.external_energy, dot(o.external_force, U))
    o.external_energy = dot(o.external_force, p.field)
    o.value = o.external_energy + o.internal_energy
    return o.value
end

# integrator hooks
function initialize!(::QuasiStaticObjectiveCache, U, p)
    p.times.time_current = zero(typeof(p.times.time_current))
    return nothing
end

function step!(solver, o::QuasiStaticObjectiveCache, U, p; verbose = true)
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p, solver.objective_cache.assembler)
    _step_begin_banner(o, p; verbose = verbose)
    solve!(solver, U, p)
    # update_field_dirichlet_bcs!(U, p.dirichlet_bcs)

    # update values at end of step
    gradient(o, U, p)
    value(o, U, p)

    # update old solution
    # update_field_dirichlet_bcs!(U, p.dirichlet_bcs)
    # copyto!(o.solution_old, U)
    copyto!(o.solution_old, p.field)

    _step_end_banner(o, p; verbose = verbose)
    return nothing
end

# logging hooks
function _step_begin_banner(::QuasiStaticObjectiveCache, p; verbose::Bool = true)
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

function _step_end_banner(o::QuasiStaticObjectiveCache, p; verbose::Bool = true)
    if verbose
        @info "External energy        = $(sum(o.external_energy))"
        @info "Internal energy        = $(sum(o.internal_energy))"
        @info "Total energy           = $(sum(o.value))"
    end
end

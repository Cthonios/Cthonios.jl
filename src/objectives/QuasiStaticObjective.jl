struct QuasiStaticObjective{
    F1 <: Function,
    F2 <: Function,
    F3 <: Function
} <: AbstractObjective{F1}
    value::F1
    gradient_u::F2
    hessian_u::F3
end

function QuasiStaticObjective()
    return QuasiStaticObjective(energy, residual, stiffness)
end

struct QuasiStaticObjectiveCache{
    # A, O, P,
    A, O,
    RT, RV <: AbstractArray{RT, 1},
    NF
} <: AbstractObjectiveCache{A, O, RT, RV}
# } <: AbstractObjectiveCache{A, O, P, RT, RV}
    assembler::A
    objective::O
    # parameters::P
    #
    external_energy::RV
    internal_energy::RV
    external_force::H1Field{RT, RV, NF}
    internal_force::H1Field{RT, RV, NF}
    # solution::H1Field{RT, RV, NF}
    solution_old::H1Field{RT, RV, NF}
    # solver helpers
    value::RV
    gradient::H1Field{RT, RV, NF}
    #
    timer::TimerOutput
end

# function QuasiStaticObjectiveCache(
#     objective::QuasiStaticObjective,
#     sim;
#     kwargs...
# )
#     # objective = QuasiStaticObjective()
#     # assembler, parameters = _setup_simulation_common(
#     #     sim, nothing; 
#     #     # return_post_processor=false,
#     #     use_condensed=true,
#     #     kwargs...
#     # )
#     assembler = _setup_assembler(
#         sim;
#         use_condensed=true,
#         kwargs...
#     )
function QuasiStaticObjectiveCache(
    assembler,
    objective::QuasiStaticObjective
)
    RT = eltype(assembler.constraint_storage)

    # backend = KA.get_backend(parameters.h1_coords)
    backend = KA.get_backend(assembler)
    external_energy = KA.zeros(backend, RT, 1)
    internal_energy = KA.zeros(backend, RT, 1)
    external_force = create_field(assembler)
    internal_force = create_field(assembler)
    # solution = create_field(assembler)
    solution_old = create_field(assembler)

    value = KA.zeros(backend, RT, 1)
    gradient = create_field(assembler)

    timer = TimerOutput()

    return QuasiStaticObjectiveCache(
        assembler, objective, #parameters,
        external_energy, internal_energy,
        external_force, internal_force,
        # solution, solution_old,
        solution_old,
        value, gradient,
        timer
    )
end

# cache hook
function setup_cache(assembler, objective)
    return QuasiStaticObjectiveCache(assembler, objective)
end

# function setup_cache(objective::O, sim; kwargs...) where O <: QuasiStaticObjective
#     return QuasiStaticObjectiveCache(objective, sim; kwargs...)
# end

# objective hooks
function gradient(o::QuasiStaticObjectiveCache, U, p)
    RT = eltype(U)
    assemble_vector!(assembler(o), o.objective.gradient_u, U, p)
    # copyto!(o.internal_force, residual(assembler(o)))
    copyto!(o.internal_force, assembler(o).residual_storage)
    fill!(o.external_force, zero(RT))
    fill!(assembler(o).residual_storage, zero(RT))
    fill!(assembler(o).residual_unknowns, zero(RT))
    assemble_vector_neumann_bc!(assembler(o), U, p)
    # copyto!(o.external_force, residual(assembler(o)))
    copyto!(o.external_force, assembler(o).residual_storage)
    o.gradient .= o.internal_force .- o.external_force

    # return FiniteElementContainers.residual(assembler(o))

    # return o.gradient.data

    if FiniteElementContainers._is_condensed(assembler(o).dof)
        FiniteElementContainers._adjust_vector_entries_for_constraints!(
          o.gradient, assembler(o).constraint_storage,
          KA.get_backend(assembler(o))
        )
        return o.gradient.data
      else
        # extract_field_unknowns!(
        #   asm.residual_unknowns, 
        #   asm.dof, 
        #   asm.residual_storage
        # )
        # return asm.residual_unknowns
        @assert false "Finish or remove me."
    end
end

function hessian(o::QuasiStaticObjectiveCache, U, p)
    assemble_stiffness!(assembler(o), o.objective.hessian_u, U, p)
    H = stiffness(assembler(o)) |> Symmetric
    return H
end

function hvp(o::QuasiStaticObjectiveCache, U, V, p)
    assemble_matrix_action!(assembler(o), o.objective.hessian_u, U, V, p)
    return FiniteElementContainers.hvp(assembler(o), V)
end

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

function value(o::QuasiStaticObjectiveCache, U, p)
    assemble_scalar!(assembler(o), o.objective.value, U, p)
    val = mapreduce(sum, sum, values(assembler(o).scalar_quadrature_storage))
    fill!(o.internal_energy, val)
    # assemble_scalar!(assembler(o), #neumann energy, U, p)

    # TODO need an actual neumann energy method in FEContainers
    fill!(o.external_energy, dot(o.external_force, U))
    o.value .= o.external_energy .+ o.internal_energy
    return sum(o.value)
end

# integrator hooks
# function initialize!(o::QuasiStaticObjectiveCache, sim)
#     # p = o.parameters
#     p = create_parameters(
#         mesh,
#         assembler, sim.physics, sim.properties; 
#         dirichlet_bcs=sim.dirichlet_bcs,
#         neumann_bcs=sim.neumann_bcs,
#         times=sim.times
#     )
#     fill!(p.times.time_current, zero(eltype(p.times.time_current)))
#     return nothing
# end
function initialize!(o::QuasiStaticObjectiveCache, U, p)
    fill!(p.times.time_current, zero(eltype(p.times.time_current)))
    return nothing
end

# function step!(o::QuasiStaticObjectiveCache, solver; verbose=true)
function step!(solver, o::QuasiStaticObjectiveCache, U, p; verbose = true)
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p)
    _step_begin_banner(o, p; verbose = verbose)
    solve!(solver, U.data, p)

    # update values at end of step
    gradient(o, U, p)
    value(o, U, p)

    # update old solution
    copyto!(o.solution_old, U)

    _step_end_banner(o, p; verbose = verbose)
    return nothing
end

# logging hooks
function _step_begin_banner(::QuasiStaticObjectiveCache, p; verbose::Bool = true)
    # p = o.parameters
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

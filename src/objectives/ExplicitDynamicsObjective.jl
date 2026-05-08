struct ExplicitDynamicsObjective{
    F1 <: Function,
    F2 <: Function,
    F3 <: Function
} <: AbstractObjective{F1}
    value::F1
    gradient_u::F2
    hessian_u::F3
end

function ExplicitDynamicsObjective()
    return ExplicitDynamicsObjective(energy, residual, stiffness)
end

mutable struct ExplicitDynamicsObjectiveCache{
    A, O,
    RT, RV <: AbstractArray{RT, 1}, NF
} <: AbstractObjectiveCache{A, O, RT, RV}
    assembler::A
    objective::O
    #
    β::RT
    γ::RT
    CFL::RT
    #
    external_energy::RV
    internal_energy::RV
    kinetic_energy::RV
    external_force::H1Field{RT, RV, NF}
    inertial_force::H1Field{RT, RV, NF}
    internal_force::H1Field{RT, RV, NF}
    solution_rate::H1Field{RT, RV, NF}
    solution_rate_rate::H1Field{RT, RV, NF}
    solution_old::H1Field{RT, RV, NF}
    solution_rate_old::H1Field{RT, RV, NF}
    #
    value::RV
    gradient::H1Field{RT, RV, NF}
    lumped_hessian::RV
    lumped_mass::RV
    #
    step_number::Int
    #
    timer::TimerOutput
end

function ExplicitDynamicsObjectiveCache(
    assembler, 
    objective::ExplicitDynamicsObjective,
    CFL, β=0.25, γ=0.5
)
    # objective = ImplicitDynamicsObjectiveNew()
    # assembler, parameters = _setup_simulation_common(
    #     sim, nothing; 
    #     return_post_processor=false,
    #     use_condensed=true
    # )

    # RT = eltype(parameters.h1_coords)
    # backend = KA.get_backend(parameters.h1_coords)
    RT = eltype(assembler.constraint_storage)
    backend = KA.get_backend(assembler)

    external_energy = KA.zeros(backend, RT, 1)
    internal_energy = KA.zeros(backend, RT, 1)
    kinetic_energy = KA.zeros(backend, RT, 1)

    external_force = create_field(assembler)
    inertial_force = create_field(assembler)
    internal_force = create_field(assembler)

    solution_rate = create_field(assembler)
    solution_rate_rate = create_field(assembler)
    
    solution_old = create_field(assembler)
    solution_rate_old = create_field(assembler)

    value = KA.zeros(backend, RT, 1)
    gradient = create_field(assembler)
    lumped_hessian = KA.zeros(backend, RT, length(gradient))
    lumped_mass = KA.zeros(backend, RT, length(gradient))

    timer = TimerOutput()

    return ExplicitDynamicsObjectiveCache(
        assembler, objective,
        β, γ, CFL,
        external_energy, internal_energy, kinetic_energy,
        external_force, inertial_force, internal_force,
        solution_rate, solution_rate_rate,
        solution_old, solution_rate_old,
        value, gradient, 
        lumped_hessian, lumped_mass,
        0, timer
    )
end

function setup_cache(assembler, objective::ExplicitDynamicsObjective, args...)
    return ExplicitDynamicsObjectiveCache(assembler, objective, args...)
end

function initialize!(
    o::ExplicitDynamicsObjectiveCache, U, p;
    lumped_mass_style = :row_sum,
    displ_ics = nothing,
    vel_ics = nothing,
    acc_ics = nothing
)
    fill!(p.times.time_current, zero(eltype(p.times.time_current)))

    M = mass_matrix(o, U, p)

    # need to lump the mass
    if lumped_mass_style == :row_sum
        n = size(M, 1)
        d = vec(sum(M, dims=2))
        o.lumped_mass .= d
    elseif lumped_mass_style == :diagonal_extraction
        # return spdiagm(0 => diag(M))
        o.lumped_mass .= diag(M)
    elseif lumped_mass_style == :scaled_diagonal
        d = diag(M)
        total_mass = sum(M)
        scale = total_mass / sum(d)
        # o.lumped_mass.= spdiagm(0 => scale .* d)
        o.lumped_mass .= scale .* d
    elseif lumped_mass_style == :block_lumping
        @assert ndim == 2 || ndim == 3
        ndofs = size(M, 1)
        @assert ndofs % ndim == 0
    
        nnodes = div(ndofs, ndim)
        Ml = zeros(Float64, ndofs)
    
        # Row-sum per dof
        rowsums = vec(sum(M, dims=2))
    
        # Accumulate nodal mass
        for a in 1:nnodes
            idx = (a-1)*ndim + 1 : a*ndim
            m_node = sum(rowsums[idx])
            Ml[idx] .= m_node / ndim
        end
    
        o.lumped_mass .= spdiagm(0 => Ml)
    else
        @assert false
    end

    # apply initial conditions 
    # TODO make this better in FEC
    if displ_ics !== nothing
        FiniteElementContainers.update_ic_values!(displ_ics, p.h1_coords)
        FiniteElementContainers.update_field_ics!(U, displ_ics)
        FiniteElementContainers.update_field_ics!(o.solution_old, displ_ics)
    end

    if vel_ics !== nothing
        FiniteElementContainers.update_ic_values!(vel_ics, p.h1_coords)
        FiniteElementContainers.update_field_ics!(o.solution_rate, vel_ics)
        FiniteElementContainers.update_field_ics!(o.solution_rate_old, vel_ics)
    end

    if acc_ics !== nothing
        FiniteElementContainers.update_ic_values!(acc_ics, p.h1_coords)
        FiniteElementContainers.update_field_ics!(o.solution_rate_rate, acc_ics)
    end

    # calculate mechanics stuff at beginning
    internal_energy = _internal_energy(o, U, p)
    kinetic_energy = 0.5 * dot(o.lumped_mass, o.solution_rate.data .* o.solution_rate.data)
    fill!(o.internal_energy, internal_energy)
    fill!(o.kinetic_energy, kinetic_energy)

    internal_force = _internal_force(o, U, p)
    o.internal_force .= internal_force

    # TODO external force
    external_force = create_field(o.assembler)
    o.external_force .= external_force
    o.inertial_force .= internal_force .- external_force
    o.solution_rate_rate.data .= o.inertial_force.data ./ o.lumped_mass
end

function _internal_energy(o::ExplicitDynamicsObjectiveCache, U, p)
    assemble_scalar!(assembler(o), o.objective.value, U, p)
    val = mapreduce(sum, +, values(assembler(o).scalar_quadrature_storage))
    fill!(o.internal_energy, val)
    # assemble_scalar!(assembler(o), #neumann energy, U, p)

    # TODO need an actual neumann energy method in FEContainers
    fill!(o.external_energy, dot(o.external_force, U))
    o.value .= o.external_energy .+ o.internal_energy
    return sum(o.value)
end

function _internal_force(o::ExplicitDynamicsObjectiveCache, U, p)
    assemble_vector!(assembler(o), o.objective.gradient_u, U, p)
    # return assembler(o).residual_storage
    o.internal_force .= assembler(o).residual_storage
end

function mass_matrix(o::ExplicitDynamicsObjectiveCache, U, p)
    # assemble_mass!(assembler(o), o.objective.hessian_u, U, p)
    # M = FiniteElementContainers.mass(assembler(o)) |> Symmetric
    # really passing displacement below instead of velocity to save
    # on having to store the solution_rate for quasi-statics
    # this method is mainly so you can construct a mass matrix
    # if you want one. Should add abstract types in the future
    assemble_mass!(assembler(o), mass, U, p)
    M = FiniteElementContainers.mass(assembler(o))
    return M #|> Symmetric
end

function step!(solver, o::ExplicitDynamicsObjectiveCache, U, p; verbose = true)
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p)
    solve!(solver, U, p)
    _step_log(o, p)
end

function _step_log(o::ExplicitDynamicsObjectiveCache, p)
    if o.step_number % 1000 == 0
        @info "Step      Time        Internal    External    Kinetic"
        @info "Number    Increment   Energy      Energy      Energy"
    end

    if o.step_number % 100 == 0
        str = @sprintf "%8d  %.4e  %.4e  %.4e  %.4e" o.step_number p.times.Δt[1] o.internal_energy[1] o.external_energy[1] o.kinetic_energy[1]
        @info str
    end

    o.step_number = o.step_number + 1
end

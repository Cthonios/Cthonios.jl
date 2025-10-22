function ExplicitDynamicsObjective()
    return QuadratureLevelObjective(energy, residual, stiffness)
end

struct ExplicitDynamicsObjectiveCache{
    A, O, P,
    RT, RV <: AbstractArray{RT, 1}, NF
} <: AbstractObjectiveCache{A, O, P, RT, RV}
    assembler::A
    objective::O
    parameters::P
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
    lumped_hessian::RV
    #
    timer::TimerOutput
end

function ExplicitDynamicsObjectiveCache(sim, CFL, β=0.25, γ=0.5)
    objective = ImplicitDynamicsObjectiveNew()
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
    lumped_hessian = KA.zeros(backend, RT, length(gradient))

    timer = TimerOutput()

    return ExplicitDynamicsObjectiveCache(
        assembler, objective, parameters,
        β, γ, CFL,
        external_energy, internal_energy, kinetic_energy,
        external_force, inertial_force, internal_force,
        predictor_solution, predictor_solution_rate,
        solution, solution_rate, solution_rate_rate,
        solution_old, solution_rate_old,
        value, gradient, lumped_hessian,
        timer
    )
end


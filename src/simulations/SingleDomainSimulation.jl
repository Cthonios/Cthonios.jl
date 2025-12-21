struct SingleDomainSimulation{
    RT <: Number, 
    RV <: AbstractArray{RT, 1},
    P1 <: NamedTuple, 
    P2 <: NamedTuple, 
    I1, 
    I2,
    I3,
    D, 
    N, 
    C
} <: AbstractSimulation
    mesh_file::String
    output_file::String
    times::TimeStepper{RV}
    physics::P1
    properties::P2
    solution_ics::I1
    solution_rate_ics::I2
    solution_rate_rate_ics::I3
    dirichlet_bcs::D
    neumann_bcs::N
    contact_pairs::C
end

function SingleDomainSimulation(
    mesh_file::String,
    output_file::String,
    times::TimeStepper,
    physics::NamedTuple,
    properties::NamedTuple;
    solution_ics::Vector{<:InitialCondition} = InitialCondition[],
    solution_rate_ics::Vector{<:InitialCondition} = InitialCondition[],
    solution_rate_rate_ics::Vector{<:InitialCondition} = InitialCondition[],
    dirichlet_bcs::Vector{<:DirichletBC} = DirichletBC[],
    neumann_bcs::Vector{<:NeumannBC} = NeumannBC[],
    contact_pairs::Vector{<:ContactPair} = ContactPair[]
)
    properties = _create_properties(physics, properties)

    return SingleDomainSimulation(
        mesh_file, output_file,
        times, physics, properties,
        solution_ics, solution_rate_ics, solution_rate_rate_ics,
        dirichlet_bcs, neumann_bcs, contact_pairs
    )
end

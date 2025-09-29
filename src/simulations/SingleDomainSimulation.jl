struct SingleDomainSimulation{M, T, P1, P2, D, N, C} <: AbstractSimulation
    mesh_file::M
    times::T
    physics::P1
    properties::P2
    dirichlet_bcs::D
    neumann_bcs::N
    contact_pairs::C
end

function SingleDomainSimulation(
    mesh_file::String,
    times::TimeStepper,
    physics::NamedTuple,
    properties::NamedTuple;
    dirichlet_bcs::Vector{<:DirichletBC} = DirichletBC[],
    neumann_bcs::Vector{<:NeumannBC} = NeumannBC[],
    contact_pairs::Vector{<:ContactPair} = ContactPair[]
)
    return SingleDomainSimulation(
        mesh_file, times, physics, properties,
        dirichlet_bcs, neumann_bcs, contact_pairs
    )
end

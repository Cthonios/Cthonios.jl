struct SingleDomainSimulation{O, U, P, M} <: AbstractSimulation
    objective::O
    u::U
    p::P
    mesh::M
    mesh_file::String
    output_file::String
end

function SingleDomainSimulation(
    objective_type,
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
    mesh = UnstructuredMesh(mesh_file)
    assembler = _setup_solid_mechanics_assembler(mesh)
    objective = objective_type(assembler)
    properties = _create_properties(physics, properties)    

    # u = create_unknowns(assembler)
    p = create_parameters(
        mesh,
        assembler, physics, properties; 
        dirichlet_bcs = dirichlet_bcs,
        ics = solution_ics,
        neumann_bcs = neumann_bcs,
        times = times
    )
    u = create_unknowns(assembler) # need to create after assembler is reset
    return SingleDomainSimulation(objective, u, p, mesh, mesh_file, output_file)
end

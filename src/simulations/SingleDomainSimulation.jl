struct SingleDomainSimulation{A, P, T} <: AbstractSimulation{A, P, T}
    assembler::A
    parameters::P
    time::T
end

function SingleDomainSimulation(
    mesh_file::String,
    times::TimeStepper,
    physics::NamedTuple,
    properties::NamedTuple;
    dirichlet_bcs::Vector{<:DirichletBC} = DirichletBC[],
    neumann_bcs::Vector{<:NeumannBC} = NeumannBC[]
)
    mesh = UnstructuredMesh(mesh_file)
    fspace = FunctionSpace(mesh, H1Field, Lagrange)

    # check consistent field names across physics
    # and setup function space elements
    funcs = map(x -> setup_function(fspace, x), values(physics))
    @assert unique(funcs) |> length == 1
    func = funcs[1]

    # setup dof manager, note all dofs are unknown at init
    if length(funcs) > 1
        dof = DofManager(func...)
    else
        dof = DofManager(func)
    end

    assembler = SparseMatrixAssembler(dof, H1Field)
    parameters = create_parameters(
        assembler, physics, properties; 
        dirichlet_bcs=dirichlet_bcs,
        neumann_bcs=neumann_bcs,
        times=times
    )
    return SingleDomainSimulation(assembler, parameters, TimerOutput())
end

using Cthonios

mesh_file = Base.source_dir() * "/cyl-tet.e"
output_file = splitext(mesh_file)[1] * "-output.exo"

times = TimeStepper(0.0, 75.0, 100)

physics = (;
    unnamed_block_1 = Cthonios.HeatConduction()
)
props = (;
    unnamed_block_1 = Dict{String, Any}(
        "density"              => 1.0,
        "heat capacity"        => 20.0,
        "thermal conductivity" => 1.0
    )
)
one_func(_, _) = 1.0
zero_func(_, _) = 0.0
dirichlet_bcs = [
    DirichletBC("temperature", zero_func; sideset_name = "bottom")
    DirichletBC("temperature", one_func; sideset_name = "top")
]

sim = SingleDomainSimulation(
    Cthonios.HeatConductionObjective,
    mesh_file, output_file,
    times, physics, props;
    dirichlet_bcs = dirichlet_bcs
)
objective = Cthonios.get_objective(sim)
u, p = Cthonios.get_unknowns_and_parameters(sim)

predictor = Cthonios.NoPredictor(sim.objective, sim.u, sim.p)
preconditioner = Cthonios.CholeskyPreconditioner(sim.objective, sim.u, sim.p)
nonlinear_solver = Cthonios.NewtonSolver(sim.objective, sim.u, sim.p)
solver = ImplicitSolver(nonlinear_solver, preconditioner, predictor)
Cthonios.run!(solver, objective, u, p, sim.mesh, output_file)

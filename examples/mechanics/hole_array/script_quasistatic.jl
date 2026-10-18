using ConstitutiveModels
using Cthonios
using FiniteElementContainers
using TimerOutputs

# file management
mesh_file = Base.source_dir() * "/mesh/hole_array.exo"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
t_end = 100.0
times = TimeStepper(0., t_end, 40)

# Physics
physics = (;
    Block1 = SolidMechanics(PlaneStrain(), Hyperelastic(NeoHookean()))
    # Block1 = SolidMechanics(PlaneStrain(), FeFv())
    # Block1 = SolidMechanics(PlaneStrain(), Hyperelastic(Hencky()))
)
# props = (;
props = Dict(
    "Block1" => Dict{String, Any}(
        "density"         => 1.,
        "Young's modulus" => 1.,
        "Poisson's ratio" => 0.48,
        # "Jm"              => 3.
        "G neq"           => 50.0,
        "relaxation time" => 50.0
    )
)

# Boundary Conditions
# func_1(x, t) = -9. * (t / t_end)
function func_1(x, t)
    if t >= t_end
        return -9.
    else
        return -9. * (t / t_end)
    end
end
func_2(x, t) = 0.0

dirichlet_bcs = [
    DirichletBC("displ_x", func_2; sideset_name = "yminus_sideset"),
    DirichletBC("displ_y", func_2; sideset_name = "yminus_sideset"),
    DirichletBC("displ_x", func_2; sideset_name = "yplus_sideset"),
    DirichletBC("displ_y", func_1; sideset_name = "yplus_sideset")
]

# Simulation
eigen_sim = SingleDomainSimulation(
    EigenObjective,
    mesh_file, output_file, 
    TimeStepper(0., 1., 1), physics, props;
    dirichlet_bcs=dirichlet_bcs
)
quasistatic_sim = SingleDomainSimulation(
    QuasiStaticObjective,
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs = dirichlet_bcs
)
# u = create_unknowns(sim.objective)
objective = Cthonios.get_objective(quasistatic_sim)
u, p = Cthonios.get_unknowns_and_parameters(quasistatic_sim)

eigen_solver = EigenSolver(eigen_sim.objective, u, p, 20)
# Cthonios.solve!(eigen_solver, eigen_sim.objective, u, p)
# vals_u = copy(eigen_sim.objective.eigen_vals)
Cthonios.run!(eigen_solver, eigen_sim.objective, u, p, eigen_sim.mesh, "eigen-undeformed.exo")

# display(vals)

timer = TimerOutput()
preconditioner = CholeskyPreconditioner(objective, u, p, timer)
predictor = TangentPredictor(objective, u, p, Val{:gmres}(), timer)
nonlinear_solver = TrustRegionSolver(objective, u, p, timer)
solver = ImplicitSolver(nonlinear_solver, preconditioner, predictor; timer = timer)
Cthonios.run!(solver, objective, u, p, quasistatic_sim.mesh, output_file)

Cthonios.run!(eigen_solver, eigen_sim.objective, u, p, eigen_sim.mesh, "eigen-deformed.exo")


# Cthonios.solve!(eigen_solver, eigen_sim.objective, u, p)
# vals_d = copy(eigen_sim.objective.eigen_vals)
# display(vals_u)
# display(vals_d)
# display(vals_u .- vals_d)

# Cthonios.run!(sim, solver)

# solver = EigenSolver(sim.objective, u, sim.p, 20; timer = timer)
# Cthonios.solve!(solver, sim.objective, u, sim.p)

#md # # Load up necessary packages
using ConstitutiveModels
using Cthonios
# using Meshes
# import Meshes: SimpleMesh, viz
# using PythonCall

#md # # File management
mesh_file = Base.source_dir() * "/mesh.g"
output_file = splitext(mesh_file)[1] * "-output.exo"

#md # # Times
times = TimeStepper(0., 1., 80)

#md # # Physics
physics = (;
    block_1 = SolidMechanics(
        ThreeDimensional(), NeoHookean()
    )
)
props = (;
    block_1 = Dict{String, Any}(
      "density"       => 1.,
      "bulk modulus"  => 10.0,
      "shear modulus" => 1.0
    )
)

#md # # Functions for BCs
func_1(x, t) = 0.0
func_2(x, t) = -0.5 * t
func_3(x, t) = @SVector [0., 0.1 * t, 0.]

#md # # Boundary conditions
dirichlet_bcs = [
  DirichletBC("displ_x", func_1; sideset_name = "sset_y_negative")
  DirichletBC("displ_y", func_1; sideset_name = "sset_y_negative")
  DirichletBC("displ_z", func_1; sideset_name = "sset_y_negative")
  DirichletBC("displ_x", func_1; sideset_name = "sset_y_positive")
  DirichletBC("displ_z", func_1; sideset_name = "sset_y_positive")
  DirichletBC("displ_y", func_2; sideset_name = "sset_y_positive")
]
neumann_bcs = NeumannBC[
  # NeumannBC("displ", func_3, "sset_y_positive")
]

#md # # Simulation setup
sim = SingleDomainSimulation(
    QuasiStaticObjective,
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs = dirichlet_bcs,
    neumann_bcs = neumann_bcs
)

#md # # objective and solution/parameters setup
objective = Cthonios.get_objective(sim)
u, p = Cthonios.get_unknowns_and_parameters(sim)

#md # # qoi extractor example
qoi = ScalarQOIExtractor(objective, energy, +, :scalar_material)

#md # # solver setup
timer = TimerOutput()
preconditioner = CholeskyPreconditioner(objective, u, p, timer)
predictor = TangentPredictor(objective, u, p, Val{:gmres}(), timer)
nonlinear_solver = TrustRegionSolver(objective, u, p, timer; max_trust_iters = 500)
solver = ImplicitSolver(nonlinear_solver, preconditioner, predictor; timer = timer)

#md # # run the simulation
Cthonios.run!(solver, objective, u, p, sim.mesh, "output.exo")

#md # # Interactive plotting
# exo = ExodusDatabase("output.e", "r")
# mesh = SimpleMesh(exo, 20)
# viz(mesh)
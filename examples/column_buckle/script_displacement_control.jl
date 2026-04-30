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
func_3(x, t) = @SVector [0., -0.025 * t, 0.]

#md # # Boundary conditions
dirichlet_bcs = [
  DirichletBC("displ_x", func_1; sideset_name = "sset_y_negative")
  DirichletBC("displ_y", func_1; sideset_name = "sset_y_negative")
  DirichletBC("displ_z", func_1; sideset_name = "sset_y_negative")
  DirichletBC("displ_x", func_1; sideset_name = "sset_y_positive")
  DirichletBC("displ_z", func_1; sideset_name = "sset_y_positive")
  DirichletBC("displ_y", func_2; sideset_name = "sset_y_positive")
]
# traction_bcs = [
#   NeumannBC("sset_y_positive", func_3)
# ]

#md # # Simulation setup
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)

#md # # objective setup
objective = Cthonios.QuasiStaticObjective()
objective_cache, U, p = Cthonios.setup_caches(objective, sim)

#md # # solver setup
solver = Cthonios.TrustRegionSolver(objective_cache, p; use_warm_start = true)
# solver = Cthonios.NewtonSolver(objective_cache)

Cthonios.run!(solver, objective_cache, U, p, sim)

#md # # Interactive plotting
# exo = ExodusDatabase("output.e", "r")
# mesh = SimpleMesh(exo, 20)
# viz(mesh)
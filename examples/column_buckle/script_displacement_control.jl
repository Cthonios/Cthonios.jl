#md # # Load up necessary packages
using ConstitutiveModels
using Cthonios
# using Meshes
# import Meshes: SimpleMesh, viz

#md # # File management
mesh_file = Base.source_dir() * "/mesh.g"
output_file = splitext(mesh_file)[1] * "-output.exo"

#md # # Times
times = TimeStepper(0., 1., 40)

#md # # Physics
physics = (;
    block_1 = SolidMechanics(
        ThreeDimensional(), NeoHookean()
    )
)
props = (;
    block_1 = Dict{String, Any}(
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
  DirichletBC("displ_x", "sset_y_negative", func_1)
  DirichletBC("displ_y", "sset_y_negative", func_1)
  DirichletBC("displ_z", "sset_y_negative", func_1)
  DirichletBC("displ_x", "sset_y_positive", func_1)
  DirichletBC("displ_z", "sset_y_positive", func_1)
  DirichletBC("displ_y", "sset_y_positive", func_2)
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
objective_cache = Cthonios.setup_cache(objective, sim)

#md # # solver setup
solver = Cthonios.TrustRegionSolverGPU(objective_cache; use_warm_start=true)

Cthonios.run!(objective_cache, solver, sim)

#md # # Interactive plotting
# exo = ExodusDatabase("output.e", "r")
# mesh = SimpleMesh(exo, 20)
# viz(mesh)
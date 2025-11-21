#md # # Load up necessary packages
using Cthonios
# using Meshes
# import Meshes: SimpleMesh, viz

#md # # File management
mesh_file = Base.source_dir() * "/mesh.g"

#md # # Functions for BCs
func_1(x, t) = 0.0
func_2(x, t) = -0.5 * t
func_3(x, t) = @SVector [0., -0.025 * t, 0.]

#md # # Boundary conditions
disp_bcs = [
  DirichletBC("nset_y_negative", [1, 2, 3], func_1)
  DirichletBC("nset_y_positive", [1, 3], func_1)
  # DirichletBC("nset_y_positive", [2], func_2)
]
traction_bcs = [
  NeumannBC("sset_y_positive", func_3)
]

#md # # Sections
sections = Section[
  Section(
    "block_1", 2,
    SolidMechanics(NeoHookean(), ThreeDimensional()),
    MaterialProperties(
      "bulk modulus" => 100.0,
      "shear modulus" => 1.0
    )
  )
]

#md # # Domain setup
domain = Domain(mesh_file, sections, disp_bcs, traction_bcs)

#md # # Objective setup
objective = Objective(domain, Cthonios.energy)

#md # # Integrator setup
integrator = QuasiStatic(0.0, 1.0, 1. / 50)

#md # # Post-processor
pp = ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y", "displ_z"])
problem = Problem(objective, integrator, TrustRegionSolver, pp; use_warm_start=false)

#md # # Finally, solve the problem
Uu, p = Cthonios.create_unknowns_and_parameters(problem)
Cthonios.solve!(problem, Uu, p)

#md # # Interactive plotting
# exo = ExodusDatabase("output.e", "r")
# mesh = SimpleMesh(exo, 20)
# viz(mesh)
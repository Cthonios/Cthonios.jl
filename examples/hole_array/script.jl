#md # # Load up necessary packages
using Cthonios
# using Meshes
import Meshes: SimpleMesh, viz

#md # # File management
mesh_file = Base.source_dir() * "/hole_array.exo"

#md # # Functions for BCs
func_1(x, t) = -5. * t#(0., -5. * t)
func_2(x, t) = 0.0
func_3(x, t) = @SVector [0., -0.025 * t]

#md # # Boundary conditions
disp_bcs = [
  DirichletBC("yminus_nodeset", [1, 2], func_2)
  DirichletBC("yplus_nodeset", [1], func_2)
  DirichletBC("yplus_nodeset", [2], func_1)
]
traction_bcs = [
]

#md # # Sections
sections = Section[
  Section(
    "Block1", 2,
    SolidMechanics(NeoHookean(), PlaneStrain()),
    MaterialProperties(
      "bulk modulus" => 0.833,
      "shear modulus" => 1.0
    )
  )
]

#md # # Domain setup
domain = Domain(mesh_file, sections, disp_bcs, traction_bcs)

#md # # Objective setup
objective = Objective(domain, Cthonios.energy)

#md # # Integrator setup
integrator = QuasiStatic(0.0, 1.0, 1. / 20)

#md # # Post-processor
pp = ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y"])
problem = Problem(objective, integrator, TrustRegionSolver, pp)

#md # # Finally, solve the problem
Uu, p = Cthonios.create_unknowns_and_parameters(problem)
Cthonios.solve!(problem, Uu, p)

#md # # Interactive plotting
exo = ExodusDatabase("output.e", "r")
mesh = SimpleMesh(exo, 20)
viz(mesh)
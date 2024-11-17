using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using LinearAlgebra
using StaticArrays
using TimerOutputs

# file management
# mesh_file = Base.source_dir() * "/hole_array.exo"
mesh_file = Base.source_dir() * "/hole_array_tri6.exo"

# functions
func_1(x, t) = -5. * t#(0., -5. * t)
func_2(x, t) = 0.0
func_3(x, t) = @SVector [0., -0.025 * t]

# bcs
disp_bcs = [
  DirichletBC("yminus_nodeset", [1, 2], func_2)
  DirichletBC("yplus_nodeset", [1], func_2)
  DirichletBC("yplus_nodeset", [2], func_1)
]
traction_bcs = [
  # NeumannBC("yplus_sideset", func_3)
]

# sections
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

# problem setup
domain = Domain(mesh_file, sections, disp_bcs, traction_bcs)
# domain = Domain(mesh_file, sections, disp_bcs)
objective = Objective(domain, Cthonios.energy)
integrator = QuasiStatic(0.0, 1.0, 1. / 80)
pp = ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y"])
problem = Problem(objective, integrator, TrustRegionSolver, pp)

# solve problem
Uu, p = Cthonios.create_unknowns_and_parameters(problem)
Cthonios.solve!(problem, Uu, p)

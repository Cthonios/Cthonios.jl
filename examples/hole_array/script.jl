using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using LinearAlgebra
using StaticArrays
using TimerOutputs

# file management
mesh_file = Base.source_dir() * "/hole_array.exo"
timer = TimerOutput()

# global setup
# times = ConstantTimeStepper(0.0, 1.0, 1. / 20)
integrator = QuasiStatic(0.0, 1.0, 1. / 20)
n_dofs = 2

# functions
func_1(x, t) = -5. * t#(0., -5. * t)
func_2(x, t) = 0.0
func_3(x, t) = @SVector [0., -0.025 * t]

# bcs
disp_bcs = [
  DirichletBC("yminus_nodeset", [1, 2], func_2)
  DirichletBC("yplus_nodeset", [1], func_2)
  # DirichletBC("yplus_nodeset", [2], func_1)
]
traction_bcs = [
  NeumannBC("yplus_sideset", func_3)
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
domain = Domain(mesh_file, sections, disp_bcs, traction_bcs, 2)
objective = Objective(
  domain, 
  Cthonios.energy, Cthonios.gradient, Cthonios.hessian, 
  Cthonios.neumann_energy, Cthonios.neumann_gradient, Cthonios.neumann_hessian, 
  timer
)
# objective = Objective(domain, Cthonios.energy, Cthonios.neumann_energy, timer)
p = ObjectiveParameters(objective, integrator)
solver = TrustRegionSolver(objective, p, timer; use_warm_start=true, preconditioner=CholeskyPreconditioner)

Uu = Cthonios.create_unknowns(solver)

# pp
pp = ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y"])

problem = QuasiStaticProblem(objective, solver, pp, timer)
Cthonios.solve!(problem, Uu, p)

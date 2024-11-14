using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using TimerOutputs

# file management
mesh_file = Base.source_dir() * "/sliding.exo"
timer = TimerOutput()

# global setup
times = ConstantTimeStepper(0.0, 1.0, 0.05)
n_dofs = 2

# functions
func_1(x, t) = 0.1 * t
func_2(x, t) = -0.1 * t
func_3(x, t) = 0.0

# bcs
disp_bcs = [
  DirichletBC("nset_left", [2], func_3),
  DirichletBC("nset_right", [2], func_3),
  DirichletBC("nset_left", [1], func_1),
  DirichletBC("nset_right", [1], func_2)
]
traction_bcs = [
]

# sections
sections = Section[
  Section(
    "unnamed_block_1", 2,
    SolidMechanics(NeoHookean(), PlaneStrain()),
    MaterialProperties(
      "bulk modulus" => 0.833,
      "shear modulus" => 0.3846
    )
  )
]

domain = Domain(mesh_file, sections, disp_bcs, traction_bcs, 2)
objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, timer)
p = ObjectiveParameters(objective, times)
solver = TrustRegionSolver(
  objective, p, timer; 
  use_warm_start=false, 
  preconditioner=CholeskyPreconditioner
)
contact = Cthonios.PenaltyContact(domain)
level_set = Cthonios.Sphere(2. / 7.)

Uu = Cthonios.create_unknowns(solver)


# pp
# pp = ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y"])

# problem = QuasiStaticProblem(objective, solver, pp, timer)
# Cthonios.solve!(problem, Uu, p)

using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using LinearAlgebra
using TimerOutputs

# file management
mesh_file = Base.source_dir() * "/window_pain_tri3.g"
timer = TimerOutput()

# global setup
times = ConstantTimeStepper(0.0, 1.0, 0.0125)
n_dofs = 2

# functions
func_1(x, t) = -0.5 * t
func_2(x, t) = 0.0

# bcs
disp_bcs = [
  DirichletBC("nset_outer_bottom", [1, 2], func_2)
  DirichletBC("nset_outer_top", [1], func_2)
  DirichletBC("nset_outer_top", [2], func_1)
]
traction_bcs = [
]

# sections
sections = TotalLagrangeSection[
  TotalLagrangeSection(
    Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
    "unnamed_block_1", 2
  )
]
domain = Domain(mesh_file, sections, disp_bcs, 2)
objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, timer)
p = ObjectiveParameters(objective, times)
# solver = NewtonSolver(objective, p, DirectSolver, timer)
solver = TrustRegionSolver(objective, p, timer; use_warm_start=false)

Uu = Cthonios.create_unknowns(solver)

# pp
pp = ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y"])

problem = QuasiStaticProblem(objective, solver, pp, timer)
Cthonios.solve!(problem, Uu, p)

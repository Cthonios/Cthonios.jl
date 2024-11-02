using ConstitutiveModels
using Cthonios
# using Enzyme
using Exodus
using FiniteElementContainers
using LinearAlgebra
using StaticArrays
using TimerOutputs

# file management
mesh_file = Base.source_dir() * "/hole_array.exo"
timer = TimerOutput()

# global setup
times = ConstantTimeStepper(0.0, 1.0, 1. / 40)
n_dofs = 2

# functions
func_1(x, t) = -5. * t#(0., -5. * t)
func_2(x, t) = 0.0
func_3(x, t) = @SVector [0., -0.0125 * t]

# bcs
disp_bcs = [
  DirichletBC("yminus_nodeset", [1, 2], func_2)
  DirichletBC("yplus_nodeset", [1], func_2)
  DirichletBC("yplus_nodeset", [2], func_1)
]
traction_bcs = [
  # NeumannBC("yplus_sideset", func_3)
]

# pre-loading
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
pp = ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y"])
objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, timer)
p = ObjectiveParameters(objective, times)
solver = TrustRegionSolver(objective, p, timer; use_warm_start=false, preconditioner=CholeskyPreconditioner)
problem = QuasiStaticProblem(objective, solver, pp, timer)
Uu = Cthonios.create_unknowns(solver)
Cthonios.solve!(problem, Uu, p)

# eigen solution

sections = Section[
  Section(
    "Block1", 2,
    SolidDynamics(NeoHookean(), PlaneStrain()),
    MaterialProperties(
      "bulk modulus" => 0.833,
      "shear modulus" => 1.0
    )
  )
]
domain = Domain(mesh_file, sections, disp_bcs, traction_bcs, 2)
pp = ExodusPostProcessor(mesh_file, "eigen_output.e", ["displ_x", "displ_y"])
objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, timer)
p = ObjectiveParameters(objective, p.t)
# solver = TrustRegionSolver(objective, p, timer; use_warm_start=false, preconditioner=CholeskyPreconditioner)
solver = EigenSolver(objective, p, timer, 10)
problem = EigenProblem(objective, solver, pp, timer)
# Uu = Cthonios.create_unknowns(solver)
Cthonios.solve!(problem, Uu, p)
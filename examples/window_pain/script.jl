using ConstitutiveModels
using Cthonios
using Enzyme
using FiniteElementContainers
using LinearSolve

# file management
mesh_file = Base.source_dir() * "/window_pain_tri3.g"

# function main()
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
    "unnamed_block_1", PlaneStrain(), NeoHookean(), 2,
    Dict{String, Any}(
      "bulk modulus" => 0.833,
      "shear modulus" => 0.3846
    )
  )
]

# domain setup
domain = Domain(
  mesh_file, times, 
  sections, disp_bcs, traction_bcs, 
  n_dofs
)

# solver
solver = NewtonSolver(domain, DirectSolver)
# solver = NewtonSolver(domain, MatrixFreeSolver)
pp = PostProcessor(
  domain, "output.e",
  ["displacement"],
  String[],
  String[]
)

# objective of problem
obj = Objective(
  Cthonios.internal_energy!,
  Cthonios.internal_force!,
  Cthonios.stiffness
)

# type of problem we want to solve
prob = Cthonios.QuasiStaticProblem(solver, obj, domain, pp)
Uu = Cthonios.solve(prob)

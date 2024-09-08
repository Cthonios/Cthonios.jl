using ConstitutiveModels
using Cthonios
using Enzyme
using FiniteElementContainers
using LinearSolve

# file management
mesh_file = Base.source_dir() * "/mesh.g"

# function main()
# global setup
times = ConstantTimeStepper(0.0, 1.0, 0.0125 / 2.)
n_dofs = 2

# functions
func_1(x, t) = 0.5 * t
func_2(x, t) = 0.0

# bcs
disp_bcs = [
  DirichletBC("nset_left", [1, 2], func_2)
  DirichletBC("nset_right", [2], func_2)
  DirichletBC("nset_right", [1], func_1)
]
traction_bcs = [
]

# sections
sections = TotalLagrangeSection[
  TotalLagrangeSection(
    "unnamed_block_1", IncompressiblePlaneStress(), Swanson(), 2,
    Dict{String, Any}(
      "bulk modulus" => 10.0,
      # "shear modulus" => 0.3846
      # "A1" => 0.264536,
      # "P1" => 3.451703,
      # "B1" => 1.89529,
      # "Q1" => 0.001,
      # "C1" => 0.856198,
      # "R1" => 0.001,
      "A1" => 0.264536,
      "alpha1" => 3.451703,
      "A2" => 0.856198,
      "alpha2" => 0.001,
      "B1" => 1.89529,
      "beta1" => 0.001,
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

try
  Uu = Cthonios.solve(prob)
catch e
  Cthonios.close(prob.post_processor)
  throw(e)
end

@show "here"
Cthonios.gradient(Reverse, prob, Uu)
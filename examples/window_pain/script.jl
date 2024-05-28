using Cthonios
using ConstitutiveModels
using Exodus

# backend for different GPUs/parallelism models
backend = NoKABackend()
common = CthoniosCommon("log.log", CthoniosBackend(backend))

# mesh
mesh = FileMesh(ExodusDatabase, "window_pain_tri3.g")
dof = DofManager(mesh)
formulation = PlaneStrain()

# functions
# funcs = ScalarFunction{2, Float64, Float64, Float64}[
#   ScalarFunction((x, t) -> 0.0, 2),
#   ScalarFunction((x, t) -> -0.5 * t, 2)
# ]
funcs = [
  (x, t) -> 0.0,
  (x, t) -> -0.5 * t
]
# materials
model = NeoHookean()
props = [
  Dict(
    "bulk modulus"  => 50.0,
    "shear modulus" => 1.0
  )
]

# sections
sections = TotalLagrangeSection[
  TotalLagrangeSection(mesh, dof, model, formulation, 1, 2)
]

# boundary conditions
disp_bcs = DisplacementBCContainer([
  DisplacementBC(mesh, 3, 1, 1),
  DisplacementBC(mesh, 3, 2, 1),
  DisplacementBC(mesh, 1, 1, 1),
  DisplacementBC(mesh, 1, 2, 2)
])
# @assert false
# time stepper
time = ConstantTimeStepper(0.0, 1.0, 0.025)

# setup a domain with all the above
domain = QuasiStaticDomain(mesh, dof, funcs, sections, props, disp_bcs, time)

# @assert false
# setup a solver
linear_solver_settings = DirectLinearSolverSettings()
linear_solver = DirectLinearSolver(domain, linear_solver_settings, backend)
solver_settings = TrustRegionSolverSettings()
solver = TrustRegionSolver(domain, linear_solver, solver_settings, backend; use_warm_start=true)

# @assert false
# setup a post-processor
pp = PostProcessor(
  "window_pain_tri3.g", "output.e",
  ["displacement"], String[], String[], 2, 2, 0, 4;
  force=true
)

# finally setup a problem and solve it
problem = ForwardProblem(domain, solver, pp)
solve!(problem, common)

# now hand off to another problem
pp = PostProcessor(
  "window_pain_tri3.g", "gradients.e",
  ["displacement", "dcoordinates"], String["dproperties"], String[], 2, 2, 0, 4;
  force=true
)
problem = EnergySensitivityProblem(problem.domain, problem.solver, pp)
solve!(problem, common)


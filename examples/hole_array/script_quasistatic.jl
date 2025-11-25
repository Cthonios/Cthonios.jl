using ConstitutiveModels
using Cthonios

# file management
mesh_file = Base.source_dir() * "/mesh/hole_array.exo"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
times = TimeStepper(0., 1., 40)

# Physics
physics = (;
    Block1 = SolidMechanics(
        PlaneStrain(), NeoHookean()
    )
)
props = (;
    Block1 = Dict{String, Any}(
        "Young's modulus" => 1.,
        "Poisson's ratio" => 0.45
    )
)

# Boundary Conditions
func_1(x, t) = -7.5 * t
func_2(x, t) = 0.0

dirichlet_bcs = [
    DirichletBC("displ_x", "yminus_sideset", func_2),
    DirichletBC("displ_y", "yminus_sideset", func_2),
    DirichletBC("displ_x", "yplus_sideset", func_2),
    DirichletBC("displ_y", "yplus_sideset", func_1)
]

# Simulation
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = Cthonios.QuasiStaticObjective()
objective_cache = Cthonios.setup_cache(objective, sim)

qoi = Cthonios.ScalarQOIExtractor(objective_cache, energy, sum, sum)

solver = Cthonios.TrustRegionSolver(objective_cache; use_warm_start=true)
# solver = Cthonios.TrustRegionSolver(objective_cache; use_warm_start=true)
# solver_type = x -> Cthonios.TrustRegionSolverGPU(x; use_warm_start=true)
Cthonios.run!(objective_cache, solver, sim) # eventually remove sim from call

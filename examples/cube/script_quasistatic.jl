using ConstitutiveModels
using Cthonios

# file management
mesh_file = Base.source_dir() * "/cube.g"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
times = TimeStepper(0., 1., 11)

# Physics
physics = (;
    cube = SolidMechanics(
        ThreeDimensional(), NeoHookean()
    )
)
props = (;
    cube = Dict{String, Any}(
        "Young's modulus" => 1.,
        "Poisson's ratio" => 0.45
    )
)

# Boundary Conditions
func_1(x, t) = 1.0 * t
func_2(x, t) = 0.0

dirichlet_bcs = [
    DirichletBC("displ_x", "ssx-", func_2),
    DirichletBC("displ_y", "ssy-", func_2),
    DirichletBC("displ_z", "ssz-", func_2),
    DirichletBC("displ_z", "ssz+", func_1)
]

# Simulation
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = Cthonios.QuasiStaticObjective()
objective_cache = Cthonios.setup_cache(objective, sim)
solver = Cthonios.TrustRegionSolver(objective_cache; use_warm_start=true)
Cthonios.run!(objective_cache, solver, sim) # eventually remove sim from call

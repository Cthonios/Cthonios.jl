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
    DirichletBC("displ_x", func_2; sideset_name = "ssx-"),
    DirichletBC("displ_y", func_2; sideset_name = "ssy-"),
    DirichletBC("displ_z", func_2; sideset_name = "ssz-"),
    DirichletBC("displ_z", func_1; sideset_name = "ssz+")
]

# Simulation
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = QuasiStaticObjective()
objective_cache, U, p = Cthonios.setup_caches(objective, sim)
solver = TrustRegionSolver(objective_cache, p; use_warm_start=true)
Cthonios.run!(solver, objective_cache, U, p, sim) # eventually remove sim from call

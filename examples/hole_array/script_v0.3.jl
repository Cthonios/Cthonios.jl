using ConstitutiveModels
using Cthonios

# Mesh
mesh_file = Base.source_dir() * "/hole_array.exo"
mesh = UnstructuredMesh(mesh_file)

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
        "bulk modulus" => 10.,
        "shear modulus" => 1.
    )
)
props = Cthonios.create_properties(physics, props)

# Boundary Conditions
func_1(x, t) = -5. * t#(0., -5. * t)
func_2(x, t) = 0.0
func_3(x, t) = @SVector [0., -0.025 * t]

dirichlet_bcs = [
    DirichletBC("displ_x", "yminus_sideset", func_2),
    DirichletBC("displ_y", "yminus_sideset", func_2),
    DirichletBC("displ_x", "yplus_sideset", func_2),
    DirichletBC("displ_y", "yplus_sideset", func_1)
]

# Simulation setup
objective = QuadratureLevelObjective(energy, residual, stiffness)
sim = SingleDomainSimulation(
    mesh_file, times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective_cache = Cthonios.QuadratureLevelObjectiveCache(objective, sim)
Uu = create_unknowns(objective_cache)
p = parameters(objective_cache)

solver = TrustRegionSolver(objective_cache, p, TimerOutput())
evolve!(objective_cache.sim_cache, solver, Uu, p)

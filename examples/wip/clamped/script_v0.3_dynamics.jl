using ConstitutiveModels
using Cthonios
using FiniteElementContainers

# Mesh
mesh_file = Base.source_dir() * "/clamped.g"

# Times
times = TimeStepper(0., 1.e-5, 11)

# Physics
physics = (;
    clamped = SolidMechanics(
        PlaneStrain(), NeoHookean()
        # PlaneStrain(), LinearElastic()
    )
)
props = (;
    clamped = Dict{String, Any}(
        "Young's modulus" => 1.,
        "Poisson's ratio" => 0.0
    )
)
props = Cthonios.create_properties(physics, props)

# Boundary Conditions
func_1(x, t) = -1. * t#(0., -5. * t)
func_2(x, t) = 0.0

dirichlet_bcs = DirichletBC[
    DirichletBC("displ_x", "yminus_sideset", func_2),
    DirichletBC("displ_y", "yminus_sideset", func_2),
    DirichletBC("displ_x", "yplus_sideset", func_2),
    DirichletBC("displ_y", "yplus_sideset", func_2)
]

# Simulation setup
objective = QuadratureLevelObjective(energy, residual, stiffness)
sim = SingleDomainSimulation(
    mesh_file, times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective_cache  = Cthonios.ImplicitDynamicsObjectiveCacheNew(sim)
solver = Cthonios.NewtonSolver(objective_cache)

# fill!(objective_cache.solution_rate, 10.)
# fill!(objective_cache.solution_rate_old, 10.)

# set initial conditions
mesh = UnstructuredMesh(mesh_file)
pp = PostProcessor(mesh, "output.e", objective_cache.assembler.dof.var)

Cthonios.run!(solver, pp)

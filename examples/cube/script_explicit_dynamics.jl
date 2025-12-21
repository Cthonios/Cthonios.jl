using ConstitutiveModels
using Cthonios
using FiniteElementContainers

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
func(_) = 1.0

vel_ics = InitialCondition[
    InitialCondition("displ_z", "cube", func)
]

# Simulation
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props
)

objective = Cthonios.ExplicitDynamicsObjective()
objective_cache, U, p = Cthonios.setup_caches(objective, sim, 0.1)

# small hack we should remove
mesh = UnstructuredMesh(mesh_file)
vel_ics = InitialConditions(mesh, objective_cache.assembler.dof, vel_ics)
# small hack above we should removes

Cthonios.initialize!(objective_cache, U, p; vel_ics = vel_ics)
solver = Cthonios.ExplicitSolver(objective_cache, p)
Cthonios.run!(solver, objective_cache, U, p, sim)

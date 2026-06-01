using ConstitutiveModels
using Cthonios
using FiniteElementContainers

# file management
mesh_file = Base.source_dir() * "/torsion.g"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
times = TimeStepper(0., 4., 11)

# Physics
physics = (;
    cube = SolidMechanics(
        ThreeDimensional(), NeoHookean()
    )
)
props = (;
    cube = Dict{String, Any}(
        "density"         => 1.e3,
        "Young's modulus" => 1.e9,
        "Poisson's ratio" => 0.25
    )
)

# Boundary Conditions
a = 8000
func_x(x) = -a * x[2] * x[3]
func_y(x) = a * x[1] * x[3]
func_z(x) = 0.0

vel_ics = InitialCondition[
    InitialCondition("displ_x", func_x; block_name = "torsion")
    InitialCondition("displ_y", func_y; block_name = "torsion")
    InitialCondition("displ_z", func_z; block_name = "torsion")
]

# Simulation
sim = SingleDomainSimulation(
    x -> Cthonios.ExplicitDynamicsObjective(x, 0.1),
    mesh_file, output_file, 
    times, physics, props
)

# small hack we should remove
mesh = UnstructuredMesh(mesh_file)
vel_ics = InitialConditions(mesh, sim.objective.assembler.dof, vel_ics)
# small hack above we should removes

Cthonios.initialize!(sim; vel_ics = vel_ics)

# @assert false
solver = Cthonios.ExplicitSolver(sim.objective, sim.p)
# Cthonios.run!(solver, objective, U, p, sim; output_exodus_every = 1000)
Cthonios.run!(sim, solver; output_exodus_every = 1000)

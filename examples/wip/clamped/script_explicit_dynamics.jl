using ConstitutiveModels
using Cthonios
using FiniteElementContainers

# Mesh
mesh_file = Base.source_dir() * "/clamped.g"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
times = TimeStepper(0., 1.e-5, 101)

# Physics
physics = (;
    clamped = SolidMechanics(
        ThreeDimensional(), LinearElastic()
    )
)
props = (;
    clamped = Dict{String, Any}(
        "density"         => 1000.0,
        "Young's modulus" => 1.0e+9,
        "Poisson's ratio" => 0.0
    )
)

# Initial conditions
a = 0.01
s = 0.02
func_1(x) = a * exp(-x[3] * x[3] / s / s / 2)

displ_ics = InitialCondition[
    InitialCondition("displ_z", "clamped", func_1)
]

# Boundary Conditions
func_2(x, t) = 0.0

dirichlet_bcs = DirichletBC[
    DirichletBC("displ_x", func_2; nodeset_name = "nsx-")
    DirichletBC("displ_x", func_2; nodeset_name = "nsx+")
    DirichletBC("displ_y", func_2; nodeset_name = "nsy-")
    DirichletBC("displ_y", func_2; nodeset_name = "nsy+")
    DirichletBC("displ_z", func_2; nodeset_name = "nsz-")
    DirichletBC("displ_z", func_2; nodeset_name = "nsz+")

]

# Simulation setup
sim = SingleDomainSimulation(
    mesh_file, output_file,
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = Cthonios.ExplicitDynamicsObjective()
objective_cache, U, p = Cthonios.setup_caches(objective, sim, 0.2)

# small hack we should remove
mesh = UnstructuredMesh(mesh_file)
displ_ics = InitialConditions(mesh, objective_cache.assembler.dof, displ_ics)
# small hack above we should removes

Cthonios.initialize!(objective_cache, U, p; displ_ics = displ_ics)
solver = Cthonios.ExplicitSolver(objective_cache, p)
Cthonios.run!(
    solver, objective_cache, U, p, sim;
    output_exodus_every = 1
)

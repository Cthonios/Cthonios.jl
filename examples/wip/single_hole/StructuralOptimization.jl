using ConstitutiveModels
using Cthonios
using FiniteElementContainers

# Setup for simulation in objective
csm_file = Base.source_dir() * "/single_hole.csm"
mesh_file = Base.source_dir() * "/single_hole_temp_temp.exo"

# Boundary conditions
func_1(x, t) = -1.5 * t#(0., -5. * t)
func_2(x, t) = 0.0

dirichlet_bcs = [
    DirichletBC("displ_x", "yminus_sideset", func_2),
    DirichletBC("displ_y", "yminus_sideset", func_2),
    DirichletBC("displ_x", "yplus_sideset", func_2),
    DirichletBC("displ_y", "yplus_sideset", func_1)
]

# Physics
physics = (;
    Block0 = Cthonios.SolidMechanics(
        PlaneStrain(), NeoHookean()
    )
)
props = (;
    Block0 = Dict{String, Any}(
        "bulk modulus" => 1.2087,
        "shear modulus" => 0.049
    )
)

props = map((x, y) -> create_properties(x, y), values(physics), values(props))
props = NamedTuple{keys(physics)}(props)

# Times
times = TimeStepper(0., 1., 20)

# objective
objective = identity

# Simulation setup
sim = SingleDomainSimulation(
    mesh_file, times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)

# optimization
opt = Cthonios.StructuralOptimization(csm_file, objective, sim)
Cthonios.optimize!(opt)

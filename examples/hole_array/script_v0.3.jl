using ConstitutiveModels
using Cthonios
using FiniteElementContainers
using StaticArrays
using TimerOutputs

# Mesh
# mesh_file = Base.source_dir() * "/hole_array.exo"
# mesh_file = "examples/window_pain_old/window_pain_tri3.g"
# mesh_file = "examples/window_pain/window_pain_temp.exo"
# mesh_file = "window.exo"
# mesh_file = "../ESP/window.exo"
mesh_file = "../ESP/hole_temp.exo"
mesh = UnstructuredMesh(mesh_file)

# Times
times = TimeStepper(0., 1., 40)

# Physics
physics = (;
    Block0 = Cthonios.SolidMechanics(
        PlaneStrain(), NeoHookean()
    )
)
props = (;
    Block0 = Dict{String, Any}(
        "bulk modulus" => 10.,
        "shear modulus" => 1.
    )
)

props = map((x, y) -> create_properties(x, y), values(physics), values(props))
props = NamedTuple{keys(physics)}(props)

# Boundary Conditions
func_1(x, t) = -2.5 * t#(0., -5. * t)
func_2(x, t) = 0.0
func_3(x, t) = @SVector [0., -0.025 * t]

dirichlet_bcs = [
    DirichletBC("displ_x", "yminus_sideset", func_2),
    DirichletBC("displ_y", "yminus_sideset", func_2),
    DirichletBC("displ_x", "yplus_sideset", func_2),
    DirichletBC("displ_y", "yplus_sideset", func_1)
    # DirichletBC("displ_x", "sset_outer_bottom", func_2),
    # DirichletBC("displ_y", "sset_outer_bottom", func_2),
    # DirichletBC("displ_x", "sset_outer_top", func_2),
    # DirichletBC("displ_y", "sset_outer_top", func_1)
]

# Solver
# solver = NewtonSolver(IterativeLinearSolver())

# Simulation setup
sim = SingleDomainSimulation(
    mesh_file, times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = Cthonios.QuadratureLevelObjective(energy, residual, stiffness)
objective_cache = Cthonios.QuadratureLevelObjectiveCache(objective, sim)
Uu = create_unknowns(objective_cache)
p = parameters(objective_cache)

# @assert false

timer = TimerOutput()
solver = Cthonios.TrustRegionSolver(objective_cache, p, timer)
displ = objective_cache.sim_cache.assembler.dof.H1_vars[1]
pp = PostProcessor(mesh, "output-hole_array.exo", displ)

for n in 1:40
    FiniteElementContainers.update_time!(p)
    FiniteElementContainers.update_bc_values!(p)
    Cthonios.solve!(solver, Uu, p)
    write_times(pp, n, FiniteElementContainers.current_time(p.times))
    write_field(pp, n, p.h1_field)
end

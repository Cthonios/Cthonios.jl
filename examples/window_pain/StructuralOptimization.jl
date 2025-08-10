using ConstitutiveModels
using Cthonios
using EngineeringSketchPadWrapper
using FiniteElementContainers

# Setup for simulation in objective
csm_file = Base.source_dir() * "/window_pain.csm"
mesh_file = Base.source_dir() * "/window_pain_temp.exo"
# mesh_file = Base.source_dir() * "/ESP_Mesh/plato/plato_CAPS.exo"

# Boundary conditions
func_1(x, t) = -0.25 * t#(0., -5. * t)
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
        "bulk modulus" => 10.,
        "shear modulus" => 1.
    )
)

props = map((x, y) -> create_properties(x, y), values(physics), values(props))
props = NamedTuple{keys(physics)}(props)

# Times
times = TimeStepper(0., 1., 40)

# Objective
objective = identity

# Simulation setup
sim = SingleDomainSimulation(
    mesh_file, times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)

# optimization
opt = Cthonios.StructuralOptimization(csm_file, objective, sim)
Cthonios.optimize!(opt)
# p = copy(opt.p)
# # p = [.15]

# results = Float64[]


# function run_gradient_check!(results, opt)
#     absteps = [1., 1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 1.e-7]
#     for step in absteps
#         g, g_fd = Cthonios.gradient_check!(opt, p, step)
#         push!(results, g_fd[1])
#     end
# end

# run_gradient_check!(results, opt)
# display(results)
# display(g_fd)
using ConstitutiveModels
using Cthonios
using EngineeringSketchPadWrapper
using FiniteElementContainers

# Setup for simulation in objective
csm_file = Base.source_dir() * "/window_pain.csm"
mesh_file = Base.source_dir() * "/window_pain_temp.exo"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Boundary conditions
func_1(x, t) = -0.25 * t
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

# Times
times = TimeStepper(0., 1., 20)

# Objective
objective = energy

# Simulation setup
sim = SingleDomainSimulation(
    mesh_file, output_file,
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)

# optimization
opt = Cthonios.StructuralOptimization(csm_file, objective, sim)
Cthonios.optimize!(opt)
# p = copy(opt.p)
# # p = [.15]

# results_ad = Float64[]
# results_fd = Float64[]


# function run_gradient_check!(results_ad, results_fd, opt)
#     # absteps = [1., 1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 1.e-7]
#     absteps = [1.e-4, 1.e-5, 1.e-6, 1.e-7, 1.e-8, 1.e-9, 5.e-9]
#     for step in absteps
#         g_ad, g_fd = Cthonios.gradient_check!(opt, p, step)
#         push!(results_ad, g_ad[1])
#         push!(results_fd, g_fd[1])
#     end
# end

# run_gradient_check!(results_ad, results_fd, opt)
# display(results_ad)
# display(results_fd)
# # display(g_fd)
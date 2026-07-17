import CthoniosDesign as CD
using ConstitutiveModels
using Cthonios

csm_file = Base.source_dir() * "/hole.csm"
opt = CD.StructuralOptimization(csm_file, :LD_MMA)

# example for how to use esp constant value
constant_p = CD.get_constant_parameters(opt)
cell_width = constant_p["cell_width"]

# simulation definition
mesh_file = Base.source_dir() * "/hole.exo"
output_file = splitext(mesh_file)[1] * "-output.exo"

times = TimeStepper(0., 1., 80)
physics = (;
    Block1 = SolidMechanics(PlaneStrain(), NeoHookean())
)
props = (;
    Block1 = Dict{String, Any}(
        "density"         => 1.,
        "Young's modulus" => 1.,
        "Poisson's ratio" => 0.48
    )
)
func_1(x, t) = -(cell_width / 3) * t
func_2(x, t) = 0.0

dirichlet_bcs = [
    DirichletBC("displ_x", func_2; sideset_name = "yminus_sideset"),
    DirichletBC("displ_y", func_2; sideset_name = "yminus_sideset"),
    DirichletBC("displ_x", func_2; sideset_name = "yplus_sideset"),
    DirichletBC("displ_y", func_1; sideset_name = "yplus_sideset")
]

sim_template = CD.SimulationTemplate(
    QuasiStaticObjective,
    mesh_file, output_file,
    times, physics, props, dirichlet_bcs
)
timer = TimerOutput()
solver = (o, u, p) -> ImplicitSolver(
    TrustRegionSolver(o, u, p, timer),
    CholeskyPreconditioner(o, u, p, timer),
    TangentPredictor(o, u, p, Val{:gmres}(), timer)
)

# TODO wrap below in    method/loop
# sim = CD.run(sim_template, solver)

qoi = ScalarQOIExtractor(energy, +, :scalar_material)
# CD.value_and_gradient!(opt, sim, qoi; last_step_only = true)

CD.optimize!(opt, sim_template, solver, qoi; maximize = true)

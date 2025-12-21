using ConstitutiveModels
using Cthonios
using Enzyme
using FiniteElementContainers
Enzyme.Compiler.VERBOSE_ERRORS[] = true

# function sim_test()
# file management
mesh_file = Base.source_dir() * "/mesh/hole_array.exo"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
# times = TimeStepper(0., 1., 4)
times = TimeStepper(0., 1., 20)

# Physics
physics = (;
    Block1 = SolidMechanics(
        PlaneStrain(), NeoHookean()
        # PlaneStrain(), LinearElastic()
        # PlaneStrain(), Gent()
        # PlaneStrain(), LinearElastoPlasticity(VonMisesYieldSurface(LinearIsotropicHardening()))
    )
)
E = 70.e3
ν = 0.3
σ_y = 200.
H = 180.
props = (;
    Block1 = Dict{String, Any}(
        "density"         => 1.,
        "Young's modulus" => 1.,
        "Poisson's ratio" => 0.48,
        # "Poisson's ratio" => 0.3,
        "Jm"              => 3.
    )
    # Block1 = Dict(
    #     "Young's modulus"           => E,
    #     "Poisson's ratio"           => ν,
    #     "yield surface"             => "VonMisesYieldSurface",
    #     "isotropic hardening model" => "LinearIsotropicHardening",
    #     "yield stress"              => σ_y,
    #     "hardening modulus"         => H
    # )
)

# Boundary Conditions
# func_1(x, t) = -0.075 * t
func_1(x, t) = -7.5 * t
func_2(x, t) = 0.0

dirichlet_bcs = [
    DirichletBC("displ_x", func_2; sideset_name = "yminus_sideset"),
    DirichletBC("displ_y", func_2; sideset_name = "yminus_sideset"),
    DirichletBC("displ_x", func_2; sideset_name = "yplus_sideset"),
    DirichletBC("displ_y", func_1; sideset_name = "yplus_sideset")
]

# Simulation
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = QuasiStaticObjective()
objective_cache, U, p = setup_caches(objective, sim)

qoi1 = QOIExtractor(
    objective_cache, helmholtz_free_energy, +,
    FiniteElementContainers.L2QuadratureField, Float64;
    is_material_qoi = true,
    reduction_2 = +
)
qoi2 = QOIExtractor(
    objective_cache, residual, identity,
    H1Field, Float64;
    is_field_qoi = true,
    reduction_2 = identity
)

solver = TrustRegionSolver(objective_cache, p; use_warm_start=true)

# solver = Cthonios.NewtonSolver(objective_cache)

Cthonios.run!(solver, objective_cache, U, p, sim)

# function design_objective!(f, qoi_storages, qois, Us, ps, params)
#     U, p = Us[end], ps[end]

#     # example for setting coords as params
#     copyto!(p.h1_coords, params)

#     # example for setting properties as params
#     # for (n, val) in enumerate(values(params))
#     #     copyto!(values(p.properties)[n], val)
#     # end
#     Cthonios.value!(f, qoi_storages[end], qois[1], U, p)
#     return nothing
# end

# obj = Cthonios.DesignObjective(design_objective!, [qoi1, qoi2], U, p)
# Cthonios.forward_problem!(obj, solver, objective_cache, U, p)

# # pick out parameters for design objective function
# coord_params = deepcopy(p.h1_coords)
# # props_params = deepcopy(p.properties)

# f = Cthonios.value(obj, coord_params)

# # f = zeros(1)
# f = create_field(objective_cache)
# Cthonios.value!(f, qoi2, U, p)
# f
# # f, dfdX = Cthonios.gradient_and_value(obj, coord_params)
# # # f, dfdp = Cthonios.gradient_and_value(obj, props_params)

# # display(f)
# # display(dfdX)
# # # # display(dfdU)
# # # # end
# # # # sim_test()
using ConstitutiveModels
using Cthonios
using FiniteElementContainers

# function sim_test()
# file management
mesh_file = Base.source_dir() * "/mesh/hole_array.exo"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
times = TimeStepper(0., 1., 40)

# Physics
physics = (;
    Block1 = SolidMechanics(
        PlaneStrain(), NeoHookean()
    )
)
props = (;
    Block1 = Dict{String, Any}(
        "Young's modulus" => 1.,
        "Poisson's ratio" => 0.45
    )
)

# Boundary Conditions
func_1(x, t) = -7.5 * t
func_2(x, t) = 0.0

dirichlet_bcs = [
    DirichletBC("displ_x", "yminus_sideset", func_2),
    DirichletBC("displ_y", "yminus_sideset", func_2),
    DirichletBC("displ_x", "yplus_sideset", func_2),
    DirichletBC("displ_y", "yplus_sideset", func_1)
]

# Simulation
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = QuasiStaticObjective()
objective_cache, U, p = setup_caches(objective, sim)

qoi = QOIExtractor(
    objective_cache, helmholtz_free_energy, +,
    FiniteElementContainers.L2QuadratureField, Float64;
    reduction_2 = +
)

solver = TrustRegionSolver(objective_cache, p; use_warm_start=true)
Cthonios.run!(solver, objective_cache, U, p, sim) # eventually remove sim from call

function design_objective!(f, qois, Us, ps, params)
    U, p = Us[end], ps[end]

    # example for setting coords as params
    copyto!(p.h1_coords, params)

    # example for setting properties as params
    # for (n, val) in enumerate(values(params))
    #     copyto!(values(p.properties)[n], val)
    # end
    Cthonios._value!(f, qois[1], U, p)
    return nothing
end

obj = Cthonios.DesignObjective(design_objective!, [qoi], U, p)

coord_params = deepcopy(p.h1_coords)
props_params = deepcopy(p.properties)
display(props_params)
push!(obj.stored_solutions, deepcopy(U))
push!(obj.stored_parameters, deepcopy(p))



# obj = Cthonios.DesignObjective(design_objective!, sens)

f = Cthonios.value(obj, coord_params)
f, dfdX = Cthonios.gradient_and_value(obj, coord_params)
# f, dfdp = Cthonios.gradient_and_value(obj, props_params)


display(f)
display(dfdX)
# # display(dfdU)
# # end
# # sim_test()
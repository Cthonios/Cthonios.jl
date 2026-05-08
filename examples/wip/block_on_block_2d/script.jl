using ConstitutiveModels
using Cthonios
using FiniteElementContainers

mesh_file = Base.source_dir() * "/mesh_tri3.exo"
output_file = splitext(mesh_file)[1] * "-output.exo"

times = TimeStepper(0., 1., 400)

physics = (;
    bottom = SolidMechanics(PlaneStrain(), NeoHookean()),
    top = SolidMechanics(PlaneStrain(), NeoHookean())
)
props = (;
    bottom = Dict{String, Any}(
        "Young's modulus" => 10.,
        "Poisson's ratio" => 0.3
    ),
    top = Dict{String, Any}(
        "Young's modulus" => 10.,
        "Poisson's ratio" => 0.3
    )
)
# props = Cthonios.create_properties(physics, props)

func_1(x, t) = 0.
func_2(x, t) = -0.5 * t

dirichlet_bcs = [
    DirichletBC("displ_x", func_1; sideset_name = "bottom"),
    DirichletBC("displ_y", func_1; sideset_name = "bottom"),
    DirichletBC("displ_x", func_1; sideset_name = "top"),
    DirichletBC("displ_y", func_2; sideset_name = "top")
]

contact_pairs = [
    ContactPair("interface_top", "interface_bottom")
]

sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective = QuasiStaticObjective()
objective_cache, U, p = setup_caches(objective, sim)
# solver = TrustRegionSolver(objective_cache, p)
# Cthonios.run!(solver, objective_cache, U, p, sim)

mesh = UnstructuredMesh(mesh_file)
dof = objective_cache.assembler.dof
# X = mesh.nodal_coords
X = p.h1_coords
contact_caches = Cthonios.create_contact_pair_caches(
    mesh, dof, contact_pairs;
    max_neighbors = 4
)

mortar = Cthonios.MortarContact(0.05, 0.1)
contact_cache = Cthonios.PenaltyContactObjectiveCache(
    mortar, contact_caches, objective_cache
)
solver = TrustRegionSolver(
    contact_cache, p; 
    max_trust_iters = 500,
    preconditioner  = Cthonios.NoPreconditioner,
    use_warm_start  = false
)
# solver = Cthonios.NewtonSolver(contact_cache)
# Cthonios.update_nearest_neighbors!(contact_caches, X, U)
# Cthonios.step!(solver, contact_cache, U, p)
# Cthonios.step!(solver, contact_cache, U, p)

# Cthonios.value(contact_cache, U, p)
# Cthonios.gradient(contact_cache, U, p)
Cthonios.run!(solver, contact_cache, U, p, sim)

# # # let's move the top block a little bit to make it penetrate
# # top_nodes = unique(getfield(mesh.element_conns, :top))
# # # @show top_nodes
# # X[2, top_nodes] .-= 0.01

# FiniteElementContainers.update_time!(p)
# FiniteElementContainers.update_bc_values!(p)
# FiniteElementContainers.update_field_dirichlet_bcs!(U, p.dirichlet_bcs)

# # Cthonios._compute_closest_distance_to_each_side(
# #     contact_caches, X, U
# # )

# function constraint_func(contact_caches, U, p)
#     Cthonios._compute_closest_distance_to_each_side(contact_caches, p.h1_coords, U)
# end

# c = constraint_func(contact_caches, U, p)
# κ = 4. * ones(length(c))
# λ = 1e-4 * abs.(κ .* c)

# # constraint_func(nothing, U, p)
# mortar = Cthonios.MortarContact(0.05, 0.1)
# penalty = Cthonios.PenaltyContact()
# al_cache = Cthonios.ConstrainedObjectiveCache(
#     objective_cache, constraint_func,
#     contact_caches
# )
# # Cthonios.gradient(al_objective_cache, U, p, λ0, κ0)
# solver = Cthonios.AugmentedLagrangeSolver(al_cache, p)
# # Cthonios.value(al_cache, U, p, λ, κ)
# # Cthonios.solve!(solver, U, p, λ, κ)

# cenergy = zeros(1)
# cvector = create_field(al_cache.objective_cache)
# Cthonios.assemble_contact_scalar!(
#     cenergy, mortar, penalty, contact_caches, X, U
# )
# Cthonios.assemble_contact_vector!(
#     cvector, mortar, penalty, contact_caches, X, U
# )

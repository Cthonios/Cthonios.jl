using ConstitutiveModels
using Cthonios
using FiniteElementContainers

mesh_file = Base.source_dir() * "/mesh_tri3.exo"

times = TimeStepper(0., 1., 20)

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
        "Young's modulus" => 0.25,
        "Poisson's ratio" => 0.3
    )
)
props = Cthonios.create_properties(physics, props)

func_1(x, t) = 0.
func_2(x, t) = -0.5 * t

dirichlet_bcs = [
    DirichletBC("displ_x", "bottom", func_1),
    DirichletBC("displ_y", "bottom", func_1),
    DirichletBC("displ_x", "top", func_1),
    DirichletBC("displ_y", "top", func_2)
]

contact_pairs = [
    ContactPair("interface_top", "interface_bottom")
]

sim = SingleDomainSimulation(
    mesh_file, times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)
objective_cache = QuasiStaticObjectiveCache(sim)
# solver = TrustRegionSolver(objective_cache, Cthonios.parameters(objective_cache), TimerOutput())
# mesh = UnstructuredMesh(mesh_file)
# pp = PostProcessor(mesh, "output.e", objective_cache.assembler.dof.var)
# Cthonios.run!(solver, pp)
# run!(solver)
# 
mesh = UnstructuredMesh(mesh_file)
X = mesh.nodal_coords
contact_enforcement = Cthonios.MortarContact(0.05, 0.1)
contact_formulation = Cthonios.PenaltyContact()
contact_caches = Cthonios.create_contact_pair_caches(mesh, contact_pairs)
# @time Cthonios.update_interactions!(contact_caches, X, U)
# @time Cthonios.update_closest_edges_and_weights!(contact_caches, X, U)
FiniteElementContainers.update_time!(objective_cache.parameters)
FiniteElementContainers.update_bc_values!(objective_cache.parameters)
U = objective_cache.solution

block_top_conns = mesh.element_conns.top.data |> unique |> sort
U[:, block_top_conns] .= -0.05
@time Cthonios.update_nearest_neighbors!(contact_caches, X, U)
# @time Cthonios.contact_potential(contact, contact_caches, X, U, 0.0005)
@time @show Cthonios.assemble_contact_scalar!(
    contact_enforcement, contact_formulation, contact_caches, X, U
)
@time Cthonios.assemble_contact_vector!(
    contact_enforcement, contact_formulation, contact_caches, X, U
)
# @time Cthonios.assemble_contact_matrix!(
#     contact_enforcement, contact_formulation, contact_caches, X, U
# )

# @show "Edge 1"
# edge_1_close = values(contact_caches)[1].interactions[:, end]
# @show X[:, values(contact_caches[1]).driver_surf.side_nodes[:, end]]

# for edge in edge_1_close
#     @show edge
#     @show X[:, values(contact_caches[1]).follow_surf.side_nodes[:, edge]]
#     @show Cthonios._compute_normal(X[:, values(contact_caches[1]).follow_surf.side_nodes[:, edge]])

# end

# add more types
struct DomainCache{
  Times, 
  Coords, Fields, BCFields,
  Props, StateOld, StateNew,
  SolverCache
}
  times::Times
  X::Coords
  U::Fields
  U_bc::BCFields
  props::Props
  state_old::StateOld
  state_new::StateNew
  solver_cache::SolverCache
  # add more fields
end

function Base.similar(cache::DomainCache)
  return DomainCache(
    similar(cache.times), similar(cache.X), similar(cache.U), 
    similar(cache.U_bc), similar(cache.props),
    similar(cache.state_old), similar(cache.state_new),
    similar(cache.solver_cache)
  )
end

struct DomainStatic{  
  Mesh,
  Dof, DBCs, NBCs,
  Sections,
  DBCDofs
}
  mesh::Mesh
  dof::Dof
  dirichlet_bcs::DBCs
  neumann_bcs::NBCs
  sections::Sections
  dirichlet_bc_dofs::DBCDofs
end

struct Domain{Static, Cache}
  static::Static
  cache::Cache
end

function Domain(
  mesh_file::String, 
  times, sections_in, 
  dirichlet_bcs_in, neumann_bcs_in,
  solver,
  n_dofs::Int
)
  mesh = FileMesh(ExodusDatabase, mesh_file)
  dof = DofManager(mesh, n_dofs)
  dirichlet_bcs = setup_dirichlet_bcs(mesh, dirichlet_bcs_in, n_dofs)
  # TODO do something about below
  neumann_bcs = NamedTuple()
  sections = setup_sections(mesh, dof, sections_in)
  Πs, props, states_old, states_new = setup_section_caches(sections_in, sections)
  static = DomainStatic(
    mesh, dof,
    dirichlet_bcs, neumann_bcs,
    sections,
    # TODO fix below
    Vector{Int64}(undef, 0) # for bc dofs
    # Vector{Float64}(undef, 0)
  )
  solver_cache = setup_solver_cache(solver, static)
  cache = DomainCache(
    times,
    coordinates(mesh),
    create_fields(dof),
    # TODO fix this type,
    Vector{Float64}(undef, 0), # for U_bc
    props,
    states_old, states_new,
    solver_cache
  )
  domain = Domain(static, cache)
  update_unknown_dofs!(domain)
  update_bcs!(domain)
  FiniteElementContainers.update_unknown_dofs!(
    domain.cache.solver_cache.assembler, 
    domain.static.dof, 
    map(x -> x.fspace, domain.static.sections), 
    domain.static.dirichlet_bc_dofs
  )
  
  return domain
end

# helpers
num_dofs_per_node(domain::Domain) = FiniteElementContainers.num_dofs_per_node(domain.static.dof)

"""
$(TYPEDSIGNATURES)
This methods assumes the bc dofs and func ids are already
properly set in the domain coming in
"""
function update_unknown_dofs!(domain::Domain)
  # update the dofs
  # TODO below line has some allocations we can minimize, but this is likely once a load step...
  all_bc_dofs = vcat(map(bc -> bc.dofs, domain.static.dirichlet_bcs)...) |> unique |> sort
  FiniteElementContainers.update_unknown_dofs!(domain.static.dof, all_bc_dofs)
  n_bc_dofs = length(all_bc_dofs)
  resize!(domain.cache.U_bc, n_bc_dofs)
  domain.cache.U_bc .= zero(eltype(domain.cache.U_bc))
  return nothing
end

function update_bcs!(domain::Domain)
  cache, static = domain.cache, domain.static
  X = cache.X
  t = cache.times.current_time
  resize!(static.dirichlet_bc_dofs, 0)
  resize!(cache.U_bc, 0)
  for bc in static.dirichlet_bcs
    func = bc.func
    for (node, dof) in zip(bc.nodes, bc.dofs)
      X_temp = @views X[:, node]
      val = func(X_temp, t)
      push!(static.dirichlet_bc_dofs, dof)
      push!(cache.U_bc, val)
    end
  end
  return nothing
end

function update_unknowns!(static, cache, Uu)
  FiniteElementContainers.update_fields!(cache.U, static.dof, Uu)
  return nothing
end

function update_fields!(static, cache, Uu)
  update_unknowns!(static, cache, Uu)
  cache.U[static.dirichlet_bc_dofs] .= cache.U_bc
  return nothing
end

function step!(domain::Domain)
  step!(domain.cache.times)
  update_bcs!(domain)
  return nothing
end

# TODO where is the best place for this
# TODO and how general should it be?
# TODO right now it's hardcoded to work for all solvers at a cost
struct SolverCache{Energy, Energies, Asm, Field1, Field2}
  Π::Energy
  Πs::Energies
  assembler::Asm
  V::Field1
  Hv::Field2
end

"""
$(TYPEDSIGNATURES)
"""
function SolverCache(dof, sections)
  assembler = StaticAssembler(dof, map(x -> x.fspace, values(sections)))
  Πs = Dict{Symbol, Any}()
  V = FiniteElementContainers.create_fields(dof)
  Hv = FiniteElementContainers.create_fields(dof)
  for (section_name, section) in pairs(sections)
    NQ, NE = FiniteElementContainers.num_q_points(section.fspace), num_elements(section.fspace)
    Πs[section_name] = ElementField{NQ, NE, Vector, Float64}(undef)
    Πs[section_name] .= zero(Float64)
  end
  return SolverCache(
    zeros(Float64, 1), ComponentArray(Πs),
    assembler, V, Hv
  )
end

"""
$(TYPEDSIGNATURES)
"""
function Base.similar(cache::SolverCache)
  Π = similar(cache.Π)
  Πs = similar(cache.Πs)
  asm = similar(cache.assembler)
  V = similar(cache.V)
  Hv = similar(cache.Hv)
  return SolverCache(Π, Πs, asm, V, Hv)
  # return SolverCache(Π, Πs, asm)
end

"""
$(TYPEDSIGNATURES)
"""
function SparseArrays.sparse!(cache::SolverCache)
  SparseArrays.sparse!(cache.assembler)
end

"""
$(TYPEDSIGNATURES)
"""
struct Domain{
  # fixed types
  Mesh, Dof, DBCs, NBCs, Sections, DBCDofs,
  # mutable types
  Times, Coords, Fields, BCFields, 
  Props,
  # Props, StateOld, StateNew, 
  SolverCache
}
  # fixed fields
  mesh::Mesh
  dof::Dof
  dirichlet_bcs::DBCs
  neumann_bcs::NBCs
  sections::Sections
  # TODO get rid of this field somehow
  dirichlet_bc_dofs::DBCDofs
  # mutable fields
  times::Times
  X::Coords
  U::Fields
  U_bc::BCFields
  props::Props
  # state_old::StateOld
  # state_new::StateNew
  solver_cache::SolverCache
end

"""
$(TYPEDSIGNATURES)
"""
function Domain(
  mesh_file::String, 
  times, sections_in, 
  dirichlet_bcs_in, neumann_bcs_in,
  # solver,
  n_dofs::Int
)
  mesh = FileMesh(ExodusDatabase, mesh_file)
  dof = DofManager(mesh, n_dofs)
  dirichlet_bcs = setup_dirichlet_bcs(mesh, dirichlet_bcs_in, n_dofs)
  # TODO do something about below
  neumann_bcs = NamedTuple()
  sections = setup_sections(mesh, dof, sections_in)
  Πs, props, states_old, states_new = setup_section_caches(sections_in, sections)
  solver_cache = SolverCache(dof, sections)
  domain = Domain(
    mesh, dof, dirichlet_bcs, neumann_bcs, sections,
    # TODO fix me
    Vector{Int64}(undef, 0), # for bc dofs
    times, coordinates(mesh), FiniteElementContainers.create_fields(dof),
    # TODO fix me
    Vector{Float64}(undef, 0), # for U_bc
    props,
    # props, states_old, states_new, 
    solver_cache
  )

  # settings dofs before returning
  update_unknown_dofs!(domain)
  update_bcs!(domain)
  FiniteElementContainers.update_unknown_dofs!(
    domain.solver_cache.assembler, 
    domain.dof, 
    map(x -> x.fspace, domain.sections), 
    domain.dirichlet_bc_dofs
  )
  return domain
end

"""
$(TYPEDSIGNATURES)
"""
function Base.similar(domain::Domain)
  return Domain(
    deepcopy(domain.mesh), deepcopy(domain.dof),
    deepcopy(domain.dirichlet_bcs), deepcopy(domain.neumann_bcs),
    deepcopy(domain.sections), deepcopy(domain.dirichlet_bc_dofs),
    #
    similar(domain.times), similar(domain.X), similar(domain.U),
    similar(domain.U_bc), 
    similar(domain.props),
    # similar(domain.state_old), similar(domain.state_new),
    similar(domain.solver_cache)
  )
end

# helpers
"""
$(TYPEDSIGNATURES)
"""
num_dofs_per_node(domain::Domain) = FiniteElementContainers.num_dofs_per_node(domain.static.dof)

"""
$(TYPEDSIGNATURES)
"""
function create_fields(domain::Domain)
  FiniteElementContainers.create_fields(domain.dof)
end

"""
$(TYPEDSIGNATURES)
"""
function create_unknowns(domain::Domain)
  FiniteElementContainers.create_unknowns(domain.dof)
end

"""
$(TYPEDSIGNATURES)
"""
function step!(domain::Domain)
  step!(domain.times)
  update_bcs!(domain)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_bcs!(domain::Domain)
  X = domain.X
  t = domain.times.current_time
  resize!(domain.dirichlet_bc_dofs, 0)
  resize!(domain.U_bc, 0)
  for bc in domain.dirichlet_bcs
    func = bc.func
    for (node, dof) in zip(bc.nodes, bc.dofs)
      X_temp = @views X[:, node]
      val = func(X_temp, t)
      push!(domain.dirichlet_bc_dofs, dof)
      push!(domain.U_bc, val)
    end
  end
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_fields!(domain, Uu)
  update_unknowns!(domain, Uu)
  domain.U[domain.dirichlet_bc_dofs] .= domain.U_bc
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_unknowns!(domain, Uu)
  FiniteElementContainers.update_fields!(domain.U, domain.dof, Uu)
  return nothing
end

"""
$(TYPEDSIGNATURES)
This methods assumes the bc dofs and func ids are already
properly set in the domain coming in
"""
function update_unknown_dofs!(domain::Domain)
  # update the dofs
  # TODO below line has some allocations we can minimize, but this is likely once a load step...
  all_bc_dofs = vcat(map(bc -> bc.dofs, domain.dirichlet_bcs)...) |> unique |> sort
  FiniteElementContainers.update_unknown_dofs!(domain.dof, all_bc_dofs)
  n_bc_dofs = length(all_bc_dofs)
  resize!(domain.U_bc, n_bc_dofs)
  domain.U_bc .= zero(eltype(domain.U_bc))
  return nothing
end

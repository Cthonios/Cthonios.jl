"""
$(TYPEDEF)
"""
abstract type AbstractDomain end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
# struct Domain{C, T, D, S, DBCs} <: AbstractDomain
struct Domain{M, D, S, DBCs} <: AbstractDomain
  # X::C
  # t::T
  mesh::M
  dof::D
  sections::S
  dirichlet_bcs::DBCs
end

"""
$(TYPEDSIGNATURES)
"""
# function Domain(mesh_file::String, times, sections_in, dbcs_in, n_dofs::Int)
function Domain(mesh_file::String, sections_in, dbcs_in, n_dofs::Int)
  mesh = FileMesh(ExodusDatabase, mesh_file)
  coords = coordinates(mesh)
  coords = NodalField{size(coords), Vector}(coords)
  dof = DofManager{n_dofs, size(coords, 2), Vector{Float64}}()
  # setup sections
  # TODO need to check all sections have compatable physics
  sections = Dict{Symbol, Any}()
  for section in sections_in
    sections[Symbol(section.block_name)] = TotalLagrangeSectionInternal(mesh, dof, section)
  end
  sections = NamedTuple(sections)
  # setup bcs
  dbcs = map(bc -> DirichletBCInternal(mesh, bc, n_dofs), dbcs_in)

  # return Domain(coords, times, dof, sections, dbcs)
  return Domain(mesh, dof, sections, dbcs)
end

"""
$(TYPEDSIGNATURES)
"""
function create_fields(domain::Domain)
  return FiniteElementContainers.create_fields(domain.dof)
end

"""
$(TYPEDSIGNATURES)
"""
function create_unknowns(domain::Domain)
  return FiniteElementContainers.create_unknowns(domain.dof)
end

"""
$(TYPEDSIGNATURES)
"""
function dirichlet_dofs(domain::Domain)
  dbcs = vcat(map(bc -> bc.dofs, domain.dirichlet_bcs)...)
  unique!(dbcs)
  sort!(dbcs)
  return dbcs
end

# """
# $(TYPEDSIGNATURES)
# """
# function step!(domain::Domain)
#   step!(domain.t)
#   return nothing
# end

"""
$(TYPEDSIGNATURES)
"""
function update_bcs!(U, domain::Domain, X, t)
  t = t.current_time
  for bc in domain.dirichlet_bcs
    for (node, dof) in zip(bc.nodes, bc.dofs)
      X_temp = @views X[:, node]
      # TODO need time updated here
      val = bc.func(X_temp, t)
      U[dof] = val
    end
  end
  return nothing
end
# function update_bcs!(U, domain::Domain)
#   t = domain.t.current_time
#   for bc in domain.dirichlet_bcs
#     for (node, dof) in zip(bc.nodes, bc.dofs)
#       X = @views domain.X[:, node]
#       # TODO need time updated here
#       val = bc.func(X, t)
#       U[dof] = val
#     end
#   end
#   return nothing
# end

"""
$(TYPEDSIGNATURES)
"""
function update_fields!(U, domain::Domain, Uu)
  FiniteElementContainers.update_fields!(U, domain.dof, Uu)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_unknown_dofs!(domain::Domain)
  dofs = dirichlet_dofs(domain)
  FiniteElementContainers.update_unknown_dofs!(domain.dof, dofs)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_unknown_dofs!(domain::Domain, asm)
  update_unknown_dofs!(domain)
  FiniteElementContainers.update_unknown_dofs!(
    asm, domain.dof, 
    map(x -> x.fspace, domain.sections),
    dirichlet_dofs(domain)
  )
  return nothing
end

"""
$(TYPEDSIGNATURES)
some FEMContainers abuse
"""
function StaticAssembler(domain::Domain)
  asm = FiniteElementContainers.StaticAssembler(
    domain.dof, map(x -> x.fspace, domain.sections)
  )
  return asm
end

# exports
export Domain
export StaticAssembler
export create_fields
export create_unknowns
export dirichlet_dofs
export update_fields!
export update_unknown_dofs!

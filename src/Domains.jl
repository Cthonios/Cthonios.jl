"""
$(TYPEDEF)
"""
abstract type AbstractDomain end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Domain{M, D, S, DBCs, DDofs} <: AbstractDomain
  mesh::M
  dof::D
  sections::S
  dirichlet_bcs::DBCs
  dirichlet_dofs::DDofs
end

"""
$(TYPEDSIGNATURES)
"""
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
  ddofs = Vector{Int}(undef, 0)
  return Domain(mesh, dof, sections, dbcs, ddofs)
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

function update_dirichlet_vals!(Ubc, domain::Domain, X, t)
  t = t.current_time
  # inefficiency below TODO
  resize!(Ubc, 0)
  for bc in domain.dirichlet_bcs
    for node in bc.nodes
      X_temp = @views X[:, node]
      val = bc.func(X_temp, t)
      push!(Ubc, val)
    end
  end
  return nothing
end

# new method below
"""
$(TYPEDSIGNATURES)
"""
function update_field_bcs!(U, domain::Domain, Ubc)
  U[domain.dirichlet_dofs] .= Ubc
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
# function update_fields!(U, domain::Domain, Uu)
function update_field_unknowns!(U, domain::Domain, Uu)
  FiniteElementContainers.update_fields!(U, domain.dof, Uu)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_unknown_dofs!(domain::Domain)
  dofs = dirichlet_dofs(domain)
  # update dof manager
  FiniteElementContainers.update_unknown_dofs!(domain.dof, dofs)
  # update dirichlet dofs
  resize!(domain.dirichlet_dofs, 0)
  for bc in domain.dirichlet_bcs
    for dof in bc.dofs
      push!(domain.dirichlet_dofs, dof)
    end
  end
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

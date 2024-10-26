"""
$(TYPEDEF)
Abstract base type for domains.
"""
abstract type AbstractDomain end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Domain type to hold information like bcs, sections, etc.
Note that the ```DofManager``` unknown dofs are not
automatically set. After setting up a ```Domain```
you will need to run, ```update_unknown_dofs!```.

```mesh``` - Handle to an open ```FileMesh```.
```dof``` - ```DofManager``` object.
```sections``` - A set of ```SectionInternal```s.
```dirichlet_bcs``` - A set of ```DirichletBCInternal```s.
```dirichlet_dofs``` - A set of dofs to apply dirichlet dofs.
  This is mainly for book-keeping purposes
"""
struct Domain{C, D, S, DBCs, DDofs, NBCs} <: AbstractDomain
  coords::C
  dof::D
  sections::S
  dirichlet_bcs::DBCs
  dirichlet_dofs::DDofs
  neumann_bcs::NBCs
end

"""
$(TYPEDSIGNATURES)
Constructor for a ```Domain``` type.
```mesh_file``` - File name of a mesh to read.
```sections_in``` - A set of ```Section```s.
```dbcs_in``` - A set of ```DirichletBC```s.
```n_dofs``` - The number of dofs in the problem.
"""
function Domain(mesh_file::String, sections_in, dbcs_in, nbcs_in, n_dofs::Int)
  mesh = FileMesh(ExodusDatabase, mesh_file)
  coords = coordinates(mesh)
  coords = NodalField{size(coords), Vector}(coords)
  dof = DofManager{n_dofs, size(coords, 2), Vector{Float64}}()
  # setup sections
  # TODO need to check all sections have compatable physics
  sections = Dict{Symbol, Any}()
  for section in sections_in
    sections[Symbol(section.block_name)] = SectionInternal(mesh, dof, section)
  end
  sections = NamedTuple(sections)
  # setup bcs
  # dbcs = map(bc -> DirichletBCInternal(mesh, bc, n_dofs), dbcs_in)
  dbcs = Dict{Symbol, Any}()
  for bc in dbcs_in
    name = bc.nset_name
    for dof in bc.dofs
      name = name * "_$dof"
    end
    dbcs[Symbol(name)] = DirichletBCInternal(mesh, bc, n_dofs)
  end
  dbcs = NamedTuple(dbcs)
  ddofs = Vector{Int}(undef, 0)

  nbcs = Dict{Symbol, Any}()
  for bc in nbcs_in
    name = bc.sset_name
    dbcs[Symbol(name)] = NeumannBCInternal(mesh, bc)
  end

  return Domain(coords, dof, sections, dbcs, ddofs, nbcs)
end

function Domain(inputs::Dict{Symbol, Any})
  mesh_file = inputs[Symbol("mesh file")]
  dbcs = inputs[Symbol("dirichlet boundary conditions")]
  dbcs = map(bc -> eval(Symbol(bc[:type]))(bc), dbcs)
  nbcs = inputs[Symbol("neumann boundary conditions")]
  nbcs = map(bc -> eval(Symbol(bc[:type]))(bc), nbcs)
  sections = inputs[:sections]
  # TODO seperate physics and sections
  sections = map(section -> eval(Symbol(section[:type]))(section), sections)
  n_dofs = map(s -> num_fields(s.physics), sections)
  @assert all(n_dofs .== n_dofs[1])
  return Domain(mesh_file, sections, dbcs, nbcs, n_dofs[1])
end

"""
$(TYPEDSIGNATURES)
Create a zero field based on ```domain.dof```.
"""
function create_fields(domain::Domain)
  return FiniteElementContainers.create_fields(domain.dof)
end

"""
$(TYPEDSIGNATURES)
Creates an unknown vector based on ```domain.dof```
"""
function create_unknowns(domain::Domain)
  return FiniteElementContainers.create_unknowns(domain.dof)
end

"""
$(TYPEDSIGNATURES)
Returns a sorted and unique vector of dirichlet dofs.
"""
function dirichlet_dofs(domain::Domain)
  dbcs = vcat(map(bc -> bc.dofs, domain.dirichlet_bcs)...)
  unique!(dbcs)
  sort!(dbcs)
  return dbcs
end

"""
$(TYPEDSIGNATURES)
Updates the values in Ubc with dirichlet boundary conditions in ```domain```.
```Ubc``` - BC values to fill.
```domain``` - Domain.
```X``` - Nodal coordinates.
```t``` - A scalar time value to use in the BC functions.
"""
function update_dirichlet_vals!(Ubc, domain::Domain, X, t)
  t = current_time(t)
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
Updates the Dirichlet BC dofs in ```U``` with the values
in ```Ubc```.
```U``` - Nodal field to update.
```domain``` - Domain.
```Ubc``` - BC values to fill.
"""
function update_field_bcs!(U, domain::Domain, Ubc)
  U[domain.dirichlet_dofs] .= Ubc
  return nothing
end

"""
$(TYPEDSIGNATURES)
Updates the unknown dofs in ```U``` with the values
in ```Uu```.
```U``` - Nodal field to update.
```domain``` - Domain.
```Uu``` - Unknown values to fill.
"""
function update_field_unknowns!(U, domain::Domain, Uu)
  FiniteElementContainers.update_fields!(U, domain.dof, Uu)
  return nothing
end

"""
$(TYPEDSIGNATURES)
Update the unknown dofs in ```domain.dof``` and
```domain.dirichlet_dofs```.
TODO maybe move ```domain.dirichlet_dofs```` to the
```DofManager``` in ```FiniteElementContainers```.
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
Update the unknown dofs in ```domain.dof```, 
```domain.dirichlet_dofs```, and ```asm```.
TODO maybe move ```domain.dirichlet_dofs```` to the
```DofManager``` in ```FiniteElementContainers```.
```domain``` - Domain object.
```asm``` - Assembly object.
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

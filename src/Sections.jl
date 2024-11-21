"""
$(TYPEDEF)
Abstract base type for sections (both for input and internal).
"""
abstract type AbstractSection end
"""
$(TYPEDSIGNATURES)
Retrieves the number of fields for the ```physics``` in ```section```
"""
num_fields(section::AbstractSection) = num_fields(section.physics)
"""
$(TYPEDSIGNATURES)
Retrieves the number of properties for the ```physics``` in ```section```
"""
num_properties(section::AbstractSection) = num_properties(section.physics)
"""
$(TYPEDSIGNATURES)
Retrieves the number of states for the ```physics``` in ```section```
"""
num_states(section::AbstractSection) = num_states(section.physics)

"""
$(TYPEDEF)
Abstract type for section inputs.
"""
abstract type AbstractSectionInput <: AbstractSection end
"""
$(TYPEDEF)
Abstract type for section internals.
```P``` corresponds to the ```physics```
```F``` corresponds to the ```formulation```
"""
abstract type AbstractSectionInternal{F, P} <: AbstractSection end

"""
$(TYPEDSIGNATURES)
Retrieves the sizes of a section in the following tuple
```(ND, NN, NP, NS)``` where
```ND``` - Number of dimensions
```NN``` - Number of nodes per element
```NP``` - Number of properties
```NS``` - Number of states
"""
function Base.size(section::AbstractSectionInternal)
  ND = FiniteElementContainers.num_dimensions(section.fspace)
  NN = FiniteElementContainers.num_nodes_per_element(section.fspace)
  NP = num_properties(section)
  NS = num_states(section)
  return ND, NN, NP, NS
end

# hack for now
function num_nodes_per_side(section::AbstractSectionInternal)
  return ReferenceFiniteElements.num_vertices(
    section.fspace.ref_fe.surface_element
  )
end

struct MaterialProperties{T}
  props::T
end

# function MaterialProperties(props::Dict{Symbol, Any})

# end

function MaterialProperties(props...)
  props = Dict{String, Any}(props)
  new_props = Dict{Symbol, Any}()
  for (k, v) in props
    new_props[Symbol(k)] = v
  end
  return MaterialProperties{typeof(new_props)}(new_props)
end

# function MaterialProperties

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Section inputs.
```physics``` - Physics object to use in this section.
```block_name``` - Name of exodus block to construct section from.
```q_order``` - Quadrature Order to use
"""
struct Section{P, M} <: AbstractSectionInput
  block_name::String
  q_order::Int
  physics::P
  props::M
end

# TODO seperate physics and sections
# TODO add properties here
function Section(inputs::Dict{Symbol, Any})
  block_name = inputs[:block]
  quadrature_order = inputs[Symbol("quadrature order")]
  physics = inputs[:physics]
  props = physics[:material][:properties]
  physics = eval(Symbol(physics[:type]))(physics)
  props = eval(Symbol(props[:type]))(props)
  return Section(block_name, quadrature_order, physics, props)
end

function element_to_block_map(mesh, sections_in)
  elem_to_block = Dict{Int, Int}()
  for sec in sections_in
    block = Block(mesh.mesh_obj, sec.block_name)
    id = block.id
    block_id_map = Exodus.read_block_id_map(mesh.mesh_obj, id)
    # display(block_id_map)
    for elem in block_id_map
      if haskey(elem_to_block, elem)
        @assert "element found in more than one block"
      end
      elem_to_block[elem] = id
    end
  end
  return elem_to_block
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Section internals.
```physics``` - Physics object to use in this section.
```block_name``` - Name of exodus block to construct section from.
```fspace``` - Function space for this element type
"""
struct SectionInternal{F, P <: AbstractPhysics, M} <: AbstractSectionInternal{F, P}
  block_id::Int
  q_order::Int
  fspace::F
  physics::P
  props::M
end

"""
$(TYPEDSIGNATURES)
Constructor for ```SectionInternal```.
```mesh``` - ```FileMesh``` object.
```dof``` - ```DofManager``` object.
```section``` - ```SectionInput``` object.
"""
function SectionInternal(mesh, dof::DofManager, section)
  block_id = mesh.mesh_obj.block_name_dict[section.block_name] |> Int
  # TODO make more efficient somehow
  elem_id_map = read_id_map(mesh.mesh_obj, ElementMap)
  elem_range = 1:Exodus.num_elements(mesh.mesh_obj.init)
  global_to_local = Dict{Int, Int}(zip(elem_id_map, elem_range))
  block_elem_id_map = read_block_id_map(
    mesh.mesh_obj, 
    mesh.mesh_obj.block_name_dict[section.block_name]
  )
  block_elem_id_map = map(x -> global_to_local[x], block_elem_id_map)
  # TODO make more efficient above
  conns = convert.(Int64, element_connectivity(mesh, section.block_name))
  conns = Connectivity{size(conns), Vector}(conns)
  elem_type = element_type(mesh, section.block_name)
  fspace = NonAllocatedFunctionSpace(dof, block_elem_id_map, conns, section.q_order, elem_type)
  props = init_properties(section.physics, section.props.props)
  return SectionInternal(block_id, section.q_order, fspace, section.physics, props)
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Section internals.
```physics``` - Physics object to use in this section.
```block_name``` - Name of exodus block to construct section from.
```fspace``` - Function space for this element type
"""
struct SurfaceSectionInternal{F, P <: AbstractPhysics, M, B} <: AbstractSectionInternal{F, P}
  block_id::Int
  q_order::Int
  fspace::F
  physics::P
  props::M
  bc::B
end

# neumann bc section
function SurfaceSectionInternal(mesh, dof::DofManager, section, bc)
  # set up book keeping stuff
  elem_ids = Vector{Int}(undef, 0)
  conns = Vector{Int}(undef, 0)
  num_nodes_per_side = Vector{Int}(undef, 0)
  for (e, element) in enumerate(bc.elements)
    # check if element in section
    if element in section.fspace.elem_id_map
      push!(elem_ids, element)
      push!(num_nodes_per_side, bc.num_nodes_per_side[e])
  
      # need to collect connectivity
      for n in 1:bc.num_nodes_per_side[e]
        k = 2 * (e - 1) + n
        push!(conns, bc.side_nodes[k])
      end
    else
      continue
    end
  end

  # check we gathered all the correct elements
  @assert all(x -> x == num_nodes_per_side[1], num_nodes_per_side)

  elem_id_map = read_id_map(mesh.mesh_obj, ElementMap)
  elem_range = 1:Exodus.num_elements(mesh.mesh_obj.init)  
  global_to_local = Dict{Int, Int}(zip(elem_id_map, elem_range))
  block_elem_id_map = read_block_id_map(mesh.mesh_obj, section.block_id)
  block_elem_id_map = map(x -> global_to_local[x], block_elem_id_map)

  # reshape conns and make it a field
  conns = Connectivity{num_nodes_per_side[1], length(elem_ids), Vector}(conns)

  # TODO fix below
  elem_type = typeof(section.fspace.ref_fe.element)
  # elem_type = eval(Symbol(String(elem_type.name.name) * "$(num_nodes_per_side[1])"))
  elem_type = eval(elem_type.name.name)
  fspace = NonAllocatedFunctionSpace(dof, block_elem_id_map, conns, section.q_order, elem_type)
  # TODO is a dummy set always the right thing?
  props = Dict{Symbol, Any}()
  return SurfaceSectionInternal(section.block_id, section.q_order, fspace, section.physics, props, bc)
end

function setup_sections(::Type{SectionInternal}, sections_in, mesh, dof, args...)
  sections = Dict{Symbol, Any}()
  for section in sections_in
    sections[Symbol(section.block_name)] = SectionInternal(mesh, dof, section)
  end
  sections = NamedTuple(sections)
  return sections
end

function setup_sections(::Type{SurfaceSectionInternal}, sections, mesh, dof, args...)
  nbcs = args[1]
  nbc_sections = Dict{Symbol, Any}()
  for (n, (section, bc)) in enumerate(Iterators.product(sections, nbcs))
    sec_name = "neumann_bc_section_$n"
    nbc_sections[Symbol(sec_name)] = SurfaceSectionInternal(mesh, dof, section, bc)
  end
  nbc_sections = NamedTuple(nbc_sections)
  return nbc_sections
end


# exports
export MaterialProperties
export Section

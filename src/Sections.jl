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

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Section internals.
```physics``` - Physics object to use in this section.
```block_name``` - Name of exodus block to construct section from.
```fspace``` - Function space for this element type
"""
struct SectionInternal{F, P <: AbstractPhysics, M} <: AbstractSectionInternal{F, P}
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
  return SectionInternal(fspace, section.physics, props)
end

# exports
export MaterialProperties
export Section

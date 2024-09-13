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
abstract type AbstractSectionInternal{P, F} <: AbstractSection end

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

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Section inputs.
```physics``` - Physics object to use in this section.
```block_name``` - Name of exodus block to construct section from.
```q_order``` - Quadrature Order to use
"""
struct Section{P} <: AbstractSectionInput
  physics::P
  block_name::String
  q_order::Int
end

# TODO seperate physics and sections
function Section(inputs::Dict{Symbol, Any})
  block_name = inputs[:block]
  quadrature_order = inputs[Symbol("quadrature order")]
  physics = inputs[:physics]
  physics = eval(Symbol(physics[:type]))(physics)
  return Section(physics, block_name, quadrature_order)
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Section internals.
```physics``` - Physics object to use in this section.
```block_name``` - Name of exodus block to construct section from.
```fspace``` - Function space for this element type
"""
struct SectionInternal{P, F} <: AbstractSectionInternal{P, F}
  physics::P
  block_name::String
  fspace::F
end

"""
$(TYPEDSIGNATURES)
Constructor for ```SectionInternal```.
```mesh``` - ```FileMesh``` object.
```dof``` - ```DofManager``` object.
```section``` - ```SectionInput``` object.
"""
function SectionInternal(mesh, dof, section)
  conns = convert.(Int64, element_connectivity(mesh, section.block_name))
  conns = Connectivity{size(conns), Vector}(conns)
  elem_type = element_type(mesh, section.block_name)
  fspace = NonAllocatedFunctionSpace(dof, conns, section.q_order, elem_type)
  return SectionInternal(section.physics, section.block_name, fspace)
end

# exports
export Section

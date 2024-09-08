"""
$(TYPEDEF)
"""
abstract type AbstractSectionInput end
"""
$(TYPEDEF)
"""
abstract type AbstractSectionInternal{P, F} end
num_fields(section::AbstractSectionInternal) = num_fields(section.physics)
num_propertie(section::AbstractSectionInternal) = num_properties(section.physics)
num_states(section::AbstractSectionInternal) = num_states(section.physics)

"""
$(TYPEDSIGNATURES)
"""
function Base.size(section::AbstractSectionInternal)
  ND = FiniteElementContainers.num_dimensions(section.fspace)
  NN = FiniteElementContainers.num_nodes_per_element(section.fspace)
  NP = num_properties(section.physics)
  NS = num_states(section.physics)
  return ND, NN, NP, NS
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Section{P} <: AbstractSectionInput
  physics::P
  block_name::String
  q_order::Int
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct SectionInternal{P, F} <: AbstractSectionInternal{P, F}
  physics::P
  block_name::String
  fspace::F
end

"""
$(TYPEDSIGNATURES)
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

# abstract type AbstractContact

abstract type AbstractContactInput end
abstract type AbstractContactInternal end

struct ContactPair{S} <: AbstractContactInput
  sset_a_name::S
  sset_b_name::S
end

struct ContactSet{E, S} <: AbstractContactInternal
  elements::E
  sides::S
  num_nodes_per_side::S
  side_nodes::S
end

function ContactSet(mesh, sset_name)
  sset = SideSet(mesh.mesh_obj, sset_name)
  elements = convert.(Int64, sset.elements)
  sides = convert.(Int64, sset.sides)
  num_nodes_per_side = convert.(Int64, sset.num_nodes_per_side)
  side_nodes = convert.(Int64, sset.side_nodes)
  return ContactSet(elements, sides, num_nodes_per_side, side_nodes)
end

struct ContactPairInternal{S1, S2} <: AbstractContactInternal
  side_a::S1
  side_b::S2
end

function ContactPairInternal(mesh, dof, section, contact_pair::ContactPair)
  sset_a = ContactSet(mesh, String(contact_pair.sset_a_name))
  sset_b = ContactSet(mesh, String(contact_pair.sset_b_name))

  sec_sset_a = SurfaceSectionInternal(mesh, dof, section, sset_a)
  sec_sset_b = SurfaceSectionInternal(mesh, dof, section, sset_b)
  # fspace_sset_a = 
  return ContactPairInternal(sec_sset_a, sec_sset_b)
end

include("LevelSets.jl")
include("PenaltyContact.jl")
include("Surfaces.jl")

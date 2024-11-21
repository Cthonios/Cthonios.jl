# TODO eventually need to get this to support multi physics
# right now it likely assumes there's a neumann condition
# on all present dofs.

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct NeumannBC{S, F} <: AbstractBCInput
  sset_name::S
  func::F
end

"""
$(TYPEDSIGNATURES)
"""
function NeumannBC(inputs::Dict{Symbol, Any})
  sideset = inputs[:sideset]
  func_expr = Meta.parse(inputs[:function])
  func = @RuntimeGeneratedFunction(func_expr)
  return NeumannBC(sideset, func)
end

function Base.show(io::IO, bc::NeumannBC)
  println(io, "NeumannBC:")
  println(io, "  Side set = $(bc.sset_name)")
  println(io, "  Function = $(bc.func)")
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct NeumannBCInternal{E, S, F} <: AbstractBCInternal
  elements::E
  sides::S
  num_nodes_per_side::S
  side_nodes::S
  func::F
end

function Base.show(io::IO, bc::NeumannBCInternal)
  println(io, "NeumannBCInternal:")
  println(io, "  Numbe of active elements = $(length(bc.elements))")
  println(io, "  Function                 = $(bc.func)")
end

"""
$(TYPEDSIGNATURES)
"""
function NeumannBCInternal(mesh, bc::NeumannBC)
  sset = SideSet(mesh.mesh_obj, bc.sset_name)
  elements = convert.(Int64, sset.elements)

  # temp_map = element_to_block_map(mesh, )

  sides = convert.(Int64, sset.sides)
  num_nodes_per_side = convert.(Int64, sset.num_nodes_per_side)
  side_nodes = convert.(Int64, sset.side_nodes)
  return NeumannBCInternal(elements, sides, num_nodes_per_side, side_nodes, bc.func)
end

# function setup_bcs(::Type{NeumannBCInternal}, mesh, nbcs_in, n_dofs)
#   nbcs = Dict{Symbol, Any}()
#   for bc in nbcs_in
#     name = bc.sset_name
#     nbcs[Symbol(name)] = NeumannBCInternal(mesh, bc)
#   end
#   nbcs = NamedTuple(nbcs)
#   return nbcs
# end

"""
$(TYPEDSIGNATURES)
"""
function setup_bcs(::Type{NeumannBCInternal}, mesh, nbcs_in, sections_in)
  elem_to_block = element_to_block_map(mesh, sections_in)
  nbcs = Dict{Symbol, Any}()

  for bc in nbcs_in
    name = bc.sset_name

    bc_temp = NeumannBCInternal(mesh, bc)
    block_ids = [elem_to_block[e] for e in bc_temp.elements]

    if !all(x -> x == block_ids[1], block_ids)
      @assert "Neumann BCs in problems with multiple element types are currently not supported"
    end

    nbcs[Symbol(name)] = NeumannBCInternal(mesh, bc)
  end
  nbcs = NamedTuple(nbcs)
  return nbcs
end
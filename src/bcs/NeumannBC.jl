# TODO eventually need to get this to support multi physics
# right now it likely assumes there's a neumann condition
# on all present dofs.

struct NeumannBC{S, F} <: AbstractBCInput
  sset_name::S
  func::F
end

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

function NeumannBCInternal(mesh, bc::NeumannBC)
  sset = SideSet(mesh.mesh_obj, bc.sset_name)
  elements = sset.elements
  sides = sset.sides
  num_nodes_per_side = sset.num_nodes_per_side
  side_nodes = sset.side_nodes
  return NeumannBCInternal(elements, sides, num_nodes_per_side, side_nodes, bc.func)
end

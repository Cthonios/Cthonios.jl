abstract type AbstractBCInput end
abstract type AbstractBCInternal end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct DirichletBC{N, D, F} <: AbstractBCInput
  nset_name::N
  dofs::D
  func::F
end

function Base.show(io::IO, bc::DirichletBC)
  println(io, "DirichletBC:")
  println(io, "  Node set = $(bc.nset_name)")
  println(io, "  Dofs     = $(bc.dofs)")
  println(io, "  Function = $(bc.func)")
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct DirichletBCInternal{N, D, F} <: AbstractBCInternal
  nodes::N
  dofs::D
  func::F
end

function Base.show(io::IO, bc::DirichletBCInternal)
  println(io, "DirichletBCInternal:")
  println(io, "  Numbe of active dofs = $(length(bc.dofs))")
  println(io, "  Function             = $(bc.func)")
end

"""
$(TYPEDSIGNATURES)
"""
function DirichletBCInternal(mesh, bc::DirichletBC, n_dofs::Int)
  for dof in bc.dofs
    if dof < 1 || dof > n_dofs
      throw(ErrorException("Bad dof number in DirichletBC."))
    end
  end

  # set up in terms of dofs, not node ids
  nset = NodeSet(mesh.mesh_obj, bc.nset_name)
  nodes = unique(nset.nodes)
  sort!(nodes)
  new_nodes = repeat(nodes, length(bc.dofs))
  new_nodes = reshape(new_nodes, length(nodes), length(bc.dofs))' |> vec

  dofs = Int64[] # Could be an issue on GPU
  for node in nodes
    for dof in bc.dofs
      push!(dofs, (node - 1) * n_dofs + dof)
    end
  end
  bc_internal = DirichletBCInternal(new_nodes, dofs, bc.func)
  return bc_internal
end

# exports
export DirichletBC

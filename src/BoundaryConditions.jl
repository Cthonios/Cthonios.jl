abstract type AbstractBC end
abstract type AbstractBCInput end

# used for inputs
struct DirichletBC{Dofs} <: AbstractBCInput
  nset_name::String
  dofs::Dofs
  func
end

function Base.show(io::IO, bc::DirichletBC)
  println(io, "DirichletBC:")
  println(io, "  Nodeset name = $(bc.nset_name)")
  println(io, "  Active dofs  = $(bc.dofs)")
  println(io, "  Function     = $(bc.func)")
end

# real driver
struct DirichletBCInternal{Nodes, Dofs, Func} <: AbstractBC
  nodes::Nodes
  dofs::Dofs
  func::Func
end

function Base.show(io::IO, bc::DirichletBCInternal)
  println(io, "DirichletBCInternal:")
  println(io, "  Numbe of active dofs = $(length(bc.dofs))")
  println(io, "  Function             = $(bc.func)")
end

# TODO add show methods for Neumann bcs

# dummy for now, no internals yet
struct NeumnanBC <: AbstractBCInput
  sset_name::String
  func
end

# TODO make Neumann bc internal do something
struct NeumnanBCInternal <: AbstractBC
  # TODO fill this out
end

# BC stuff
"""
Add a single BC. Note that update_unknown_dofs will have to be called after this
"""
function setup_dirichlet_bc(mesh, bc::DirichletBC, n_dofs)
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

"""
Add a list of bcs. update_unknown_dofs is always called at the end
even if there are no dirichlet bcs
"""
function setup_dirichlet_bcs(
  mesh, bcs_in::Vector{T}, n_dofs
) where T <: DirichletBC
  bcs = Dict{Symbol, Any}()
  for (n, bc) in enumerate(bcs_in)
    bc_name = Symbol("dirichlet_bc_$n")
    bcs[bc_name] = setup_dirichlet_bc(mesh, bc, n_dofs)
  end
  return NamedTuple(bcs)
end

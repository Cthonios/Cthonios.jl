"""
$(TYPEDEF)
$(TYPEDFIELDS)
Base Dirichlet boundary condition type used 
for inputs either from a script or input file. 
```nset_name``` corresponds to the name of the 
node set in the exodus file, ```dofs``` correspond
to the indexed fields this bc is to be applied to,
and ```func``` is the function to apply to the fields.
"""
struct DirichletBC{N, D, F} <: AbstractBCInput
  nset_name::N
  dofs::D
  func::F
end

function DirichletBC(inputs::Dict{Symbol, Any})
  nodeset = inputs[:nodeset]
  dofs = inputs[:dofs]
  func_expr = Meta.parse(inputs[:function])
  func = @RuntimeGeneratedFunction(func_expr)
  return DirichletBC(nodeset, dofs, func)
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
Base internal Dirichlet boundary condition used for
internal purposes.
```nodes``` corresponds to the node ids this BC is to be
applied to,
```dofs``` is the set of degrees of freedoms this bc 
is to be applied to,
and ```func``` is the function to apply to a field on the dofs.
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
Constructor for internal Dirichlet boundary condition.
```mesh``` is the exodus mesh to read node sets from,
```bc``` is the ```DirichletBC``` input, and ```n_dofs``` is
the total number of fields in the problem.
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
  new_nodes = reshape(new_nodes, length(nodes), length(bc.dofs))' |> vec |> collect 

  dofs = Int64[] # Could be an issue on GPU
  for node in nodes
    for dof in bc.dofs
      push!(dofs, (node - 1) * n_dofs + dof)
    end
  end
  bc_internal = DirichletBCInternal(new_nodes, dofs, bc.func)
  return bc_internal
end

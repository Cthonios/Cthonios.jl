"""
"""
struct DisplacementBC{Nodes, Dofs} <: AbstractDirichletBC{Nodes, Dofs}
  nodes::Nodes
  dofs::Dofs
  func_id::Int64
  # func::ScalarFunction{2, Float64, Float64, Float64} # Change if this works
end

"""
"""
function DisplacementBC(inputs::D, mesh_file::FileMesh, num_dofs::Int, nset_id::Int, dof::Int) where D <: Dict{Symbol, Any}
  func_id = inputs[:function][:func_id]
  nset_nodes = convert.(Int64, nodeset(mesh_file, nset_id))
  sort!(nset_nodes)
  dofs = similar(nset_nodes)
  for n in axes(dofs, 1)
    dofs[n] = (nset_nodes[n] - 1) * num_dofs + dof
  end
  # func_ids = keys(funcs)
  @info "  Displacement BC"
  @info "    Nodeset ID   = $nset_id"
  @info "    Dof          = $dof"
  @info "    Function     = $(inputs[:function][:expression])"

  return DisplacementBC{typeof(nset_nodes), typeof(dofs)}(nset_nodes, dofs, func_id)
end

"""
Function to help in script inputs

Note this assumes the number of dofs is equal to the number of dimensions
int he mesh file. This will break for thermomechanical problems

TODO change DisplacementBC to hold a unique FunctionWrapper
TODO also add a method that uses a nset name
"""
function DisplacementBC(mesh_file::FileMesh, nset_id::Int, dof::Int, func_id::Int)
  num_dofs = FiniteElementContainers.num_dimensions(mesh_file)
  nset_nodes = convert.(Int64, nodeset(mesh_file, nset_id))
  sort!(nset_nodes)
  dofs = similar(nset_nodes)
  for n in axes(dofs, 1)
    dofs[n] = (nset_nodes[n] - 1) * num_dofs + dof
  end
  return DisplacementBC{typeof(nset_nodes), typeof(dofs)}(nset_nodes, dofs, func_id)
end

# gives IR barf
# function DisplacementBC(mesh_file::FileMesh, nset_id::Int, dof::Int, func)
#   num_dofs = FiniteElementContainers.num_dimensions(mesh_file)
#   nset_nodes = convert.(Int64, nodeset(mesh_file, nset_id))
#   sort!(nset_nodes)
#   dofs = similar(nset_nodes)
#   for n in axes(dofs, 1)
#     dofs[n] = (nset_nodes[n] - 1) * num_dofs + dof
#   end
#   return DisplacementBC{typeof(nset_nodes), typeof(dofs)}(nset_nodes, dofs, func)
# end

function setup_displacement_bcs(
  inputs::D, mesh_file::FileMesh, num_dofs::Int
) where D <: Vector{Dict{Symbol, Any}}
# )::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}} where D <: Vector{Dict{Symbol, Any}}

  new_section("Displacement Boundary Conditions")
  @info "Reading in and setting up displacement BCs"

  # collect displacement bcs from input file
  bcs = DisplacementBC[]
  for bc in inputs
    for nset_id in bc[Symbol("nodeset ids")]
      for dof in bc[Symbol("dofs")]
        push!(bcs, DisplacementBC(bc, mesh_file, num_dofs, nset_id, dof))
      end
    end
  end

  # return collect_displacement_bc_indices(bcs)
  return DisplacementBCContainer(bcs)
end

function collect_displacement_bc_indices(
  bcs
)
  bc_nodes = mapreduce(x -> x.nodes, vcat, bcs)
  bc_dofs = mapreduce(x -> x.dofs, vcat, bcs)
  bc_func_ids = mapreduce(x -> fill(x.func_id, length(x.dofs)), vcat, bcs)

  # filter out repeat dofs, this is a first come first serve
  unique_bc_dof_ids = unique(i -> bc_dofs[i], eachindex(bc_dofs))

  unique!(bc_dofs)
  # TODO maybe add sorting to the node ids for faster access?
  bc_nodes = bc_nodes[unique_bc_dof_ids]
  bc_func_ids = bc_func_ids[unique_bc_dof_ids]

  # sort them
  sort_bc_nodes = sortperm(bc_nodes)

  bc_nodes = bc_nodes[sort_bc_nodes]
  bc_dofs = bc_dofs[sort_bc_nodes]
  bc_func_ids = bc_func_ids[sort_bc_nodes]

  return bc_nodes, bc_dofs, bc_func_ids
end

"""
"""
struct DisplacementBCContainer{
  BCNodes,
  BCDofs,
  BCFuncIDs
}
  bc_nodes::BCNodes
  bc_dofs::BCDofs
  bc_func_ids::BCFuncIDs
end

"""
"""
function DisplacementBCContainer(bcs)
  bc_nodes, bc_dofs, bc_func_ids = collect_displacement_bc_indices(bcs)
  return DisplacementBCContainer(bc_nodes, bc_dofs, bc_func_ids)
end

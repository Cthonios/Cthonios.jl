struct DisplacementBC{Nodes, Dofs} <: AbstractDirichletBC{Nodes, Dofs}
  nodes::Nodes
  dofs::Dofs
  func_id::Int64
end

function DisplacementBC(inputs::D, mesh_file::FileMesh, num_dofs::Int) where D <: Dict{Symbol, Any}
  nset_id = inputs[Symbol("nodeset id")]
  dof = inputs[:dof]
  # func_inputs = inputs[:function][]
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
  # @info "    Nodeset Name = $nset_name"
  @info "    Dof          = $dof"
  @info "    Function     = $(inputs[:function][:expression])"
  # func = Symbol(bc["function"])
  # func_id = findall(x -> x == func, func_ids)[1]
  # return setup_function(func_inputs)
  # func = setup_function(func_inputs)

  return DisplacementBC{typeof(nset_nodes), typeof(dofs)}(nset_nodes, dofs, func_id)
end

function setup_displacement_bcs(
  inputs::D, mesh_file::FileMesh, num_dofs::Int
)::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}} where D <: Vector{Dict{Symbol, Any}}

  new_section("Displacement Boundary Conditions")
  @info "Reading in and setting up displacement BCs"
  bcs = tuple(map(x -> DisplacementBC(x, mesh_file, num_dofs), inputs)...)

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

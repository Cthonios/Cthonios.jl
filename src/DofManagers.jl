struct DofManager{B <: AbstractArray{Bool, 2}, V <: AbstractArray{<:Integer}}
  is_unknown::B
  unknown_indices::V
end

"""
By default sets it up to have no BCs

User must update BCs later
"""
function DofManager(coords::M) where {M <: AbstractMatrix}
  # hard-coded to solid mechanics
  n_dofs, n_nodes = size(coords)

  # get total element count by summing up over all blocks
  # n_els = map(x -> size(x.conn, 2), blocks) |> sum

  # indexing arrays
  is_unknown      = BitArray(1 for _ = 1:n_dofs, _ = 1:n_nodes)
  ids             = reshape(1:length(is_unknown), n_dofs, n_nodes)
  unknown_indices = ids[is_unknown]

  return DofManager(is_unknown, unknown_indices)
end

"""
"""
create_fields(d::DofManager, Rtype::Type{<:AbstractFloat} = Float64) = zeros(Rtype, size(d.is_unknown))

"""
"""
create_unknowns(d::DofManager, Rtype::Type{<:AbstractFloat} = Float64) = zeros(Rtype, sum(d.is_unknown))

"""
"""
dof_ids(d::DofManager) = reshape(1:length(d.is_unknown), size(d.is_unknown))

"""
"""
Base.size(d::DofManager) = size(d.is_unknown)

# NOTE this is only for a static case
# Time 0. is fed into here
function update_displacement_bcs!(
  dof::DofManager, U::M1, bcs::V, coords::M2, t::F
) where {M1 <: AbstractMatrix, M2 <: AbstractMatrix, V <: AbstractVector{<:EssentialBC}, F <: AbstractFloat}

  dof.is_unknown .= 1
  # for bc in bcs
  #   for node in bc.nodes
  #     @time dof.is_unknown[bc.dof, node] = 0
  #     @time U[bc.dof, node] = @views bc.func(coords[:, node], t)
  #   end
  # end
  for bc in bcs
    dof.is_unknown[bc.dof, bc.nodes] .= 0
    U[bc.dof, bc.nodes] = @views bc.func.(eachcol(coords[:, bc.nodes]), (t,))
  end
end

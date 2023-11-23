# TODO eventually move this code to FiniteElementContainers
 
# connectivity containers
abstract type AbstractConnectivity{N, NDof, NxNDof} end

function setup_static_connectivity(block::B, n_nodes::Int, NDof::Int; type::Type{<:Integer} = Int64) where B <: Block
  N   = size(block.conn, 1)
  ids = reshape(1:NDof * n_nodes, NDof, n_nodes)
  connectivity = reinterpret(SVector{N, Int64}, vec(convert.(type, block.conn)))
  dof_conn = ids[:, block.conn]
  dof_conn = reinterpret(SVector{N * NDof, Int64}, vec(dof_conn))
  conn     = StructArray(connectivity)
  dof_conn = StructArray(dof_conn)
  return conn, dof_conn
end

struct Connectivity{N, NDof, NxNDof, A1, A2} <: AbstractConnectivity{N, NDof, NxNDof}
  conn::A1
  dof_conn::A2
end

# come up with a cleaner way to do below  TODO
function Connectivity(block::B, n_nodes::Int, n_dofs::Int, type::Type = SVector) where {B <: Block}
  N = size(block.conn, 1)
  NxNDof = N * n_dofs
  if type <: SVector
    conn, dof_conn = setup_static_connectivity(block, n_nodes, n_dofs)
  else
    error("Not supported")
  end

  return Connectivity{N, n_dofs, NxNDof, typeof(conn), typeof(dof_conn)}(conn, dof_conn)
end

connectivity(conn::Connectivity, e::Int) = conn.conn[e]
dof_connectivity(conn::Connectivity, e::Int) = conn.dof_conn[e]

###################################################################################

# These containers are meant to capture
# all of the interpolation done in fem calculations

abstract type AbstractFunctionSpace{N, D, NDof, NxNDof, RefFE, Conn} end
num_nodes_per_element(::F) where {
  F <: AbstractFunctionSpace{N, D, NDof, NxNDof, RefFE, Conn} where {N, D, NDof, NxNDof, RefFE, Conn}
} = N
num_dimensions(::F) where {
  F <: AbstractFunctionSpace{N, D, NDof, NxNDof, RefFE, Conn} where {N, D, NDof, NxNDof, RefFE, Conn}
} = D
reference_element(fspace::F) where F <: AbstractFunctionSpace = fspace.ref_fe
connectivity(fspace::F) where F <: AbstractFunctionSpace = connectivity(fspace.conn)
connectivity(fspace::F, e::Int) where F <: AbstractFunctionSpace = connectivity(fspace.conn, e)
dof_connectivity(fspace::F) where F <: AbstractFunctionSpace = dof_connectivity(fspace.conn)
dof_connectivity(fspace::F, e::Int) where F <: AbstractFunctionSpace = dof_connectivity(fspace.conn, e)

number_of_q_points(fspace::F) where F <: AbstractFunctionSpace = length(fspace.ref_fe.interpolants)

function map_shape_function_gradients(X::V1, ∇N_ξ::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}
  J     = X * ∇N_ξ
  J_inv = inv(J)
  ∇N_Xs = (J_inv * ∇N_ξ')'
  return ∇N_Xs, det(J)
end

function setup_preallocated_interpolants!(Xs, ∇N_Xs, JxWs, el_coords, re)
  for e in axes(Xs, 2)
    for q in axes(Xs, 1)
      N = ReferenceFiniteElements.shape_function_values(re, q)
      Xs[q, e] = el_coords[e] * N
      ∇N_X, detJ  = map_shape_function_gradients(el_coords[e], ReferenceFiniteElements.shape_function_gradients(re, q))
      ∇N_Xs[q, e] = ∇N_X
      JxWs[q, e]  = quadrature_weight(re, q) * detJ
    end
  end
end

struct PreAllocatedFunctionSpace{
  N, D, NDof, NxNDof, 
  RefFE <: ReferenceFE, Conn <: AbstractConnectivity,
  S1, S2, S3
} <: AbstractFunctionSpace{N, D, NDof, NxNDof, RefFE, Conn}

  ref_fe::RefFE
  conn::Conn
  Xs::S1
  ∇N_Xs::S2
  JxWs::S3
end

function PreAllocatedFunctionSpace(coords::M,
  block::B,
  n_dofs::Int,
  q_degree::Int,
  type::Type = Float64
) where {M <: AbstractMatrix, B <: Block}

  ref_fe = ReferenceFE(block, q_degree)
  conn   = Connectivity(block, size(coords, 2), n_dofs)

  N      = size(block.conn, 1)
  D      = size(coords, 1)
  NDof   = n_dofs
  NxNDof = N * NDof
  Q      = length(ref_fe.interpolants)
  E      = size(block.conn, 2)

  el_coords = reinterpret(SMatrix{D, N, type, D * N}, vec(@views coords[:, block.conn]))
  Xs    = StructArray{SVector{D, type}}(undef, Q, E)
  ∇N_Xs = StructArray{SMatrix{N, D, type, N * D}}(undef, Q, E)
  JxWs  = Matrix{type}(undef, Q, E)
  setup_preallocated_interpolants!(Xs, ∇N_Xs, JxWs, el_coords, ref_fe)
  return PreAllocatedFunctionSpace{N, D, NDof, NxNDof, typeof(ref_fe), typeof(conn), typeof(Xs), typeof(∇N_Xs), typeof(JxWs)}(
    ref_fe, conn, Xs, ∇N_Xs, JxWs
  )
end

Base.show(io::IO, ::PreAllocatedFunctionSpace) = println(io, "PreAllocatedFunctionSpace")

number_of_elements(fspace::PreAllocatedFunctionSpace) = size(fspace.JxWs, 2)

quadrature_point(fspace::PreAllocatedFunctionSpace, q::Int, e::Int) = fspace.Xs[q, e]
quadrature_points(fspace::PreAllocatedFunctionSpace)                = fspace.Xs

shape_function_values(fspace::PreAllocatedFunctionSpace, q::Int, ::Int) = 
ReferenceFiniteElements.shape_function_values(reference_element(fspace), q)

# TODO find a way to reduce allocations below. 
# maybe using lazy or fill arrays?
# function shape_function_values(fspace::PreAllocatedFunctionSpace)
#   # repeat(ReferenceFiniteElements.shape_function_values(reference_element(fspace)), 
#   #       number_of_elements(fspace)) |> StructArray
#   # temp = repate(ReferenceFE)
#   Ns = ReferenceFiniteElements.shape_function_values(reference_element(fspace))
#   temp = repeat(Ns, )
# end

shape_function_gradients(fspace::PreAllocatedFunctionSpace, q::Int, e::Int) = fspace.∇N_Xs[q, e]
shape_function_gradients(fspace::PreAllocatedFunctionSpace)                 = fspace.∇N_Xs

# """
# function space that does not pre-allocate
# all of the shape function values
# but rather calculates them on the fly
# from the reference finite element
# """
# struct NonAllocatedFunctionSpace{
#   N, D, RefFE <: ReferenceFE, Conn <: AbstractConnectivity
# } <: AbstractFunctionSpace{N, D, RefFE, Conn}
#   ref_fe::RefFE
#   conn::Conn
# end

# function NonAllocatedFunctionSpace(
#   coords::M,
#   block::B,
#   n_dofs::Int,
#   q_degree::Int
# ) where {M <: AbstractMatrix, B <: Block}

#   num_nodes, num_elem = size(block.conn)
#   num_dims            = size(coords, 1)

#   # conn, dof_conn = setup_connectivities(coords, block)
#   ref_fe = ReferenceFE(block, q_degree) 
#   conn   = Connectivity(block, size(coords, 2), n_dofs)
#   return NonAllocatedFunctionSpace{num_nodes, num_dims, typeof(ref_fe), typeof(conn)}(
#     ref_fe, conn
#   )
# end

# shape_function_values(fspace::F, q::Int, ::Int) where F <: NonAllocatedFunctionSpace =
# ReferenceFiniteElements.shape_function_values(fspace.ref_fe, q)

# function shape_function_gradients(fspace::F, q::Int, e::Int) where F <: NonAllocatedFunctionSpace

# end

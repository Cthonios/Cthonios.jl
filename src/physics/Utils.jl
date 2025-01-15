# helper methods
function reshape_coordinates(::AbstractPhysics, cell, X_el)
  ND = size(cell.∇N_ξ, 2)
  NN = length(X_el) ÷ ND
  return SMatrix{ND, NN, eltype(X_el), ND * NN}(X_el...)
end

"""
comes in as (N_nodes x n_fields) vector
returns as n_fields x n_nodes matrix
"""
function reshape_field(physics::AbstractPhysics, cell, U_el)
  NF = num_fields(physics)
  # NN = length(cell.N)
  NN = length(U_el) ÷ NF
  # TODO is this division necessary?
  return SMatrix{NF, NN, eltype(U_el), NF * NN}(U_el...)
end

function interpolate_field_values(physics::AbstractPhysics, cell, U_el)
  U_el = reshape_field(physics, cell, U_el)
  return U_el * cell.N
end

function interpolate_field_gradients(physics::AbstractPhysics, cell, U_el)
  U_el = reshape_field(physics, cell, U_el)
  return U_el * cell.∇N_X
end

function interpolate_field_values_and_gradients(physics::AbstractPhysics, cell, U_el)
  U_el = reshape_field(physics, cell, U_el)
  return U_el * cell.N, U_el * cell.∇N_X
end

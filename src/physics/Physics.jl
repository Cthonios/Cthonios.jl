"""
$(TYPEDEF)
```NF``` - Number of fields in this physics\n
```NP``` - Number of properties in this physics\n
```NS``` - Number of states in this physics
"""
abstract type AbstractPhysics{NF, NP, NS} end

"""
$(TYPEDSIGNATURES)
"""
num_fields(::AbstractPhysics{NF, NP, NS}) where {NF, NP, NS} = NF
"""
$(TYPEDSIGNATURES)
"""
num_properties(::AbstractPhysics{NF, NP, NS}) where {NF, NP, NS} = NP
"""
$(TYPEDSIGNATURES)
"""
num_states(::AbstractPhysics{NF, NP, NS}) where {NF, NP, NS} = NS



# helper methods
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

# implementations
include("Poisson.jl")
include("SolidMechanics.jl")

# top level energy, etc.
# these methods below essentially convert things from FEM
# lingo to stuff that can be passed to kernels that look
# more like the math
function energy(physics::AbstractPhysics, cell, u_el, times, state, props)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q, ∇u_q = interpolate_field_values_and_gradients(physics, cell, u_el)
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  val = energy(physics, u_q, ∇u_q, X_q, props)
  return JxW * val
end

function gradient(physics::AbstractPhysics, cell, u_el, times, state, props)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q, ∇u_q = interpolate_field_values_and_gradients(physics, cell, u_el)
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  G = discrete_gradient(physics.formulation, ∇N_X)
  val = gradient(physics, u_q, ∇u_q, N, G, X_q, props)
  return JxW * val
end

function hessian(physics::AbstractPhysics, cell, u_el, times, state, props)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q, ∇u_q = interpolate_field_values_and_gradients(physics, cell, u_el)
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  G = discrete_gradient(physics.formulation, ∇N_X)
  val = hessian(physics, u_q, ∇u_q, N, G, X_q, props)
  return JxW * val
end

# TODO how to overload so we can fall back to AD methods when
# e.g. gradient and hessian are not defined out of laziness
# or difficulty (same thing really)

# Neumann bc methods
function energy(physics::AbstractPhysics, bc::NeumannBCInternal, cell, u_el, times)
  @unpack X_q, N, N_reduced, JxW, n = cell
  # u_q = interpolate_field_values(physics, cell, u_el)
  u_el = reshape_field(physics, cell, u_el)
  # return U_el * cell.N
  u_q = u_el * cell.N_reduced
  t_q = bc.func(X_q, times.current_time)
  return -JxW * dot(u_q, t_q)
end

function gradient(::AbstractPhysics, bc::NeumannBCInternal, cell, u_el, times)
  @unpack X_q, N, JxW, n = cell
  t_q = bc.func(X_q, current_time(times))
  f_q = t_q * N'
  return -JxW * f_q[:]
end

function hessian(physics::AbstractPhysics, bc::NeumannBCInternal, cell, u_el, times)
  N = length(cell.N) * num_fields(physics)
  return zeros(SMatrix{N, N, eltype(u_el), N * N})
end

# exports
export Poisson
export SolidMechanics

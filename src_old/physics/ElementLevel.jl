# top level energy, etc.
# these methods below essentially convert things from FEM
# lingo to stuff that can be passed to kernels that look
# more like the math
function energy(physics::AbstractPhysics, cell, u_el, X_el, state, props, times)
  # cell = MappedInterpolants(cell, )
  X_el = reshape_coordinates(physics, cell, X_el)
  cell = MappedInterpolants(cell, X_el)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q, ∇u_q = interpolate_field_values_and_gradients(physics, cell, u_el)
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  val = energy(physics, u_q, ∇u_q, X_q, times, state, props)
  return JxW * val
end

function gradient(physics::AbstractPhysics, cell, u_el, X_el, state, props, times)
  X_el = reshape_coordinates(physics, cell, X_el)
  cell = MappedInterpolants(cell, X_el)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q, ∇u_q = interpolate_field_values_and_gradients(physics, cell, u_el)
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  G = discrete_gradient(physics.formulation, ∇N_X)
  val = gradient(physics, u_q, ∇u_q, N, G, X_q, times, state, props)
  return JxW * val
end

# function gradient(physics::AbstractPhysics, cell, u_el, times, state, props)
#   return ForwardDiff.gradient(z -> energy(physics, cell, z, times, state, props), u_el)
# end

function hessian(physics::AbstractPhysics, cell, u_el, X_el, state, props, times)
  X_el = reshape_coordinates(physics, cell, X_el)
  cell = MappedInterpolants(cell, X_el)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q, ∇u_q = interpolate_field_values_and_gradients(physics, cell, u_el)
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  G = discrete_gradient(physics.formulation, ∇N_X)
  val = hessian(physics, u_q, ∇u_q, N, G, X_q, times, state, props)
  return JxW * val
end

function mass_matrix(physics::AbstractPhysics, cell, u_el, X_el, state, props, times)
  X_el = reshape_coordinates(physics, cell, X_el)
  cell = MappedInterpolants(cell, X_el)
  @unpack X_q, N, ∇N_X, JxW = cell
  u_q, ∇u_q = interpolate_field_values_and_gradients(physics, cell, u_el)
  ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
  N = FiniteElementContainers.discrete_values(physics.formulation, N)
  G = discrete_gradient(physics.formulation, ∇N_X)
  val = mass_matrix(physics, u_q, ∇u_q, N, G, X_q, times, state, props)
  return JxW * val
end

# TODO need to fix and patch up a lot of below

# TODO how to overload so we can fall back to AD methods when
# e.g. gradient and hessian are not defined out of laziness
# or difficulty (same thing really)

# Neumann bc methods
function neumann_energy(physics::AbstractPhysics, cell, u_el, X_el, times, ::NeumannBCInternal, bc_val, ref_fe, q, face)
  X_el = reshape_coordinates(physics, cell, X_el)
  cell = MappedSurfaceInterpolants(ref_fe, X_el, q, face)
  @unpack X_q, N, N_reduced, JxW, n = cell
  u_el = reshape_field(physics, cell, u_el)
  u_q = u_el * cell.N_reduced
  # u_q = u_el * cell.N
  # t_q = bc.func(X_q, current_time(times)) # allocations here
  t_q = bc_val
  return -JxW * dot(u_q, t_q)
end

function neumann_gradient(physics::AbstractPhysics, cell, u_el, X_el, times, ::NeumannBCInternal, bc_val, ref_fe, q, face)
  # X_el = reshape_coordinates(physics, cell, X_el)
  # cell = MappedSurfaceInterpolants(cell, )
  X_el = reshape_coordinates(physics, cell, X_el)
  cell = MappedSurfaceInterpolants(ref_fe, X_el, q, face)
  @unpack X_q, N, N_reduced, JxW, n = cell
  # t_q = bc.func(X_q, current_time(times)) # allocations here
  t_q = bc_val
  # f_q = t_q * N'
  f_q = t_q * N_reduced'
  return -JxW * f_q[:]
end

# TODO this is incorrect for the current interface below
# Needs to return a SMatrix of size N_reduced
function neumann_hessian(physics::AbstractPhysics, cell, u_el, X_el, times, ::NeumannBCInternal, bc_val, ref_fe, q, face)
  N = length(cell.N) * num_fields(physics)
  return zeros(SMatrix{N, N, eltype(u_el), N * N})
end

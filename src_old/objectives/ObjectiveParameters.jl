"""
$(TYPEDEF)
$(TYPEDFIELDS)
Type for objective function parameters for design parameters
such as coordinates, time, bc values, properties
state variables, and some scratch arrays.
"""
struct ObjectiveParameters{U1, T, B, N, S, P, U3, V1, U4, Q, IC} <: AbstractObjectiveParameters
# struct ObjectiveParameters{U1, T, B, N, S, P, U2, U4, Q} <: AbstractObjectiveParameters
  # design parameters
  X::U1
  t::T # TODO eventually move out since now it's an integrator
  Ubc::B # dirichlet bc design parameters
  nbc::N # neumann bc design parameters
  state_old::S
  state_new::S
  props::P
  # scratch arrays
  # U::U2
  grad_scratch::U3
  grad_vec_scratch::V1
  hvp_scratch::U4
  q_vals_scratch::Q
  # caches
  integrator_cache::IC
end

"""
$(TYPEDSIGNATURES)
Constructor for a ```ObjectiveParameters``` type.
```o``` - Objective function object.
```times``` - Times object.
"""
# TODO change below to have times renamed to integrator
function ObjectiveParameters(o::AbstractObjective, times)
  X = copy(o.domain.coords)
  # U = create_fields(o.domain)
  # boundary conditions
  Ubc = Vector{eltype(X)}(undef, 0)
  nbc = Vector{SVector{size(X, 1), eltype(X)}}(undef, 0)
  update_dirichlet_vals!(Ubc, o.domain, X, times)
  update_neumann_vals!(nbc, o.domain, X, times)
  # properties
  props = map(sec -> sec.props, o.domain.sections)
  props = ComponentArray(props)
  # scratch arrays
  grad_scratch = create_fields(o.domain)
  grad_vec_scratch = create_unknowns(o.domain)
  hvp_scratch = create_fields(o.domain)
  # need a scratch array for calculating q values on gpus
  # TODO move somewhere else
  q_vals_scratch = Dict{Symbol, Any}()
  state_old = Dict{Symbol, Any}()
  state_new = Dict{Symbol, Any}()
  for (name, sec) in pairs(o.domain.sections)
    NQ = FiniteElementContainers.num_q_points(sec.fspace)
    NE = FiniteElementContainers.num_elements(sec.fspace)
    q_vals_scratch[name] = zeros(eltype(X), NQ, NE)
    state_old[name] = repeat(init_state(sec.physics), outer=(1, NQ, NE))
    state_new[name] = repeat(init_state(sec.physics), outer=(1, NQ, NE))
  end
  q_vals_scratch = ComponentArray(q_vals_scratch)
  state_old = ComponentArray(state_old)
  state_new = ComponentArray(state_new)

  # in development stuff
  int_cache = integrator_cache(o, times)

  params = ObjectiveParameters(
    X, times, Ubc, nbc, state_old, state_new, props,
    grad_scratch, grad_vec_scratch, hvp_scratch, q_vals_scratch,
    int_cache
    # U, hvp_scratch, q_vals_scratch
  )
  return params
end

# bindings to integrator methods
current_solution(p::ObjectiveParameters) = current_solution(p.integrator_cache)

# TODO can probably replace this with make_zero! from enzyme
function zero_parameters!(p::ObjectiveParameters)
  p.X .= zero(eltype(p.X))
  # p.t .= 
  if length(p.Ubc) > 0
    p.Ubc .= zero(eltype(p.Ubc))
  end
  if length(p.nbc) > 0
    for n in axes(p.nbc, 1)
      p.nbc[n] = zero(typeof(p.nbc[n]))
    end
  end
  p.state_old .= zero(eltype(p.state_old))
  p.state_new .= zero(eltype(p.state_new))
  p.props .= zero(eltype(p.props))
  # p.U .= zero(eltype(p.U))
  U = current_solution(p)
  U .= zero(eltype(U))
  p.hvp_scratch .= zero(eltype(p.hvp_scratch))
  p.q_vals_scratch .= zero(eltype(p.q_vals_scratch))
  return nothing
end

# methods for updating design parameters
"""
$(TYPEDSIGNATURES)
"""
function update_dirichlet_vals!(p::ObjectiveParameters, o::AbstractObjective)
  update_dirichlet_vals!(p.Ubc, o.domain, p.X, p.t)
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function update_neumann_vals!(p::ObjectiveParameters, o::AbstractObjective)
  update_neumann_vals!(p.nbc, o.domain, p.X, p.t)
  return nothing
end

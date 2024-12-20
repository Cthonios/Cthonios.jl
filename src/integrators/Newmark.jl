struct Newmark{I, T, U} <: AbstractTimeIntegrator
  # time
  start_time::T
  end_time::T
  current_time::U
  current_time_step::I
  Δt::U
  # newmark parameters
  β::T
  γ::T
end

"""
$(TYPEDSIGNATURES)
Method to construct a ```Newmark```.
```start_time``` - the initial time value.
```end_time``` - time to end the simulation.
```Δt``` - the time step to use for all time steps.
```β``` - Newmark β parameter
```γ``` - Newmark γ parameter
"""
function Newmark(start_time::T, end_time::T, Δt::T, β::T, γ::T) where T <: Number
  return Newmark(start_time, end_time, [start_time], [1], [Δt], β, γ)
end

function Newmark(inputs::Dict{Symbol, Any})
  start_time = inputs[Symbol("start time")]
  end_time = inputs[Symbol("end time")]
  Δt = inputs[Symbol("time increment")]
  β = inputs[Symbol("beta")]
  γ = inputs[Symbol("gamma")]
  return Newmark(start_time, end_time, Δt, β, γ)
end

"""
$(TYPEDSIGNATURES)
"""
function integration_step_header(times::Newmark)
  @info "$(repeat('=', 96))"
  @info "= Load step    $(times.current_time_step[1])"
  @info "= Old Time     $(times.current_time[1])"
  @info "= New Time     $(times.current_time[1] + times.Δt[1])"
  @info "= End Time     $(times.end_time[1])"
  @info "= % Completete $(100.0 * (times.current_time[1] + times.Δt[1]) / times.end_time[1])"
  @info "$(repeat('=', 96))"
end

"""
$(TYPEDSIGNATURES)
Method to increment ```time.current_time``` by ```Δt```.
"""
function step!(time::Newmark) 
  temp = current_time(time) + time_step(time)
  time.current_time[1] = temp
  time.current_time_step[1] += 1
  return nothing
end

# cache
struct NewmarkCache{Us}
  U::Us
  V::Us
  A::Us
  U_old::Us
  V_old::Us
  A_old::Us
end

function NewmarkCache(objective)
  U = create_fields(objective.domain)
  V = create_fields(objective.domain)
  A = create_fields(objective.domain)
  U_old = create_fields(objective.domain)
  V_old = create_fields(objective.domain)
  A_old = create_fields(objective.domain)
  return NewmarkCache(U, V, A, U_old, V_old, A_old)
end

integrator_cache(objective, ::Newmark) = NewmarkCache(objective)

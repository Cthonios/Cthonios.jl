struct QuasiStatic{I, T, U} <: AbstractQuasiStaticTimeIntegrator
  start_time::T
  end_time::T
  current_time::U
  current_time_step::I
  Δt::U
end

current_time(t::QuasiStatic) = t.current_time[1]
end_time(t::QuasiStatic) = t.end_time
start_time(t::QuasiStatic) = t.start_time
time_step(t::QuasiStatic) = t.Δt[1]

function QuasiStatic(inputs::Dict{Symbol, Any})
  start_time = inputs[Symbol("start time")]
  end_time = inputs[Symbol("end time")]
  Δt = inputs[Symbol("time increment")]
  return QuasiStatic(start_time, end_time, Δt)
end

"""
$(TYPEDSIGNATURES)
"""
function integration_step_header(times::QuasiStatic)
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
function step!(time::QuasiStatic) 
  temp = current_time(time) + time_step(time)
  time.current_time[1] = temp
  time.current_time_step[1] += 1
  return nothing
end

"""
$(TYPEDSIGNATURES)
Method to construct a ```QuasiStatic```.
```start_time``` - the initial time value.
```end_time``` - time to end the simulation.
```Δt``` - the time step to use for all time steps.
"""
function QuasiStatic(start_time::T, end_time::T, Δt::T) where T <: Number
  return QuasiStatic(start_time, end_time, [start_time], [1], [Δt])
end

function Base.show(io::IO, time::QuasiStatic)
  println(io, "QuasiStatic:")
  println(io, "  Start time = $(time.start_time)")
  println(io, "  End time   = $(time.end_time)")
  println(io, "  Time step  = $(time.Δt)")
end

# function scratch_arrays(objective, ::QuasiStatic)
#   U = create_fields(objective.domain)
# end
# scratch arrays necessary for the integrator to work
# these go in the parameters object p
struct QuasiStaticCache{T <: NodalField}
  U::T
end

function QuasiStaticCache(objective)
  U = create_fields(objective.domain)
  return QuasiStaticCache{typeof(U)}(U)
end

current_solution(cache::QuasiStaticCache) = cache.U

integrator_cache(objective, ::QuasiStatic) = QuasiStaticCache(objective)

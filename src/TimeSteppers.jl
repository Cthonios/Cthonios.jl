"""
$(TYPEDEF)
Abstract type for time steppers. 
Expects a type ```T``` to correspond to the type
for time values such as ```Float64``` or a time
with a unit from for example ```Unitful```.
"""
abstract type AbstractTimeStepper{T} end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
Time stepper type with a constant time step.
```start_time``` - the initial time value.
```end_time``` - time to end the simulation.
```current_time``` - the current time stored by this time stepper.
```current_time_step``` - the current index of the time step.
```Δt``` - the time step to use for all time steps.
"""
mutable struct ConstantTimeStepper{T <: Number} <: AbstractTimeStepper{T}
  start_time::T
  end_time::T
  current_time::T
  current_time_step::Int
  Δt::T
end

"""
$(TYPEDSIGNATURES)
Method to increment ```time.current_time``` by ```Δt```.
"""
function step!(time::ConstantTimeStepper) 
  temp = time.current_time + time.Δt
  time.current_time = temp
  time.current_time_step += 1
  return nothing
end

"""
$(TYPEDSIGNATURES)
Method to construct a ```ConstantTimeStepper```.
```start_time``` - the initial time value.
```end_time``` - time to end the simulation.
```Δt``` - the time step to use for all time steps.
"""
function ConstantTimeStepper(start_time::T, end_time::T, Δt::T) where T <: Number
  return ConstantTimeStepper(start_time, end_time, start_time, 1, Δt)
end

function Base.show(io::IO, time::ConstantTimeStepper)
  println(io, "ConstantTimeStepper:")
  println(io, "  Start time = $(time.start_time)")
  println(io, "  End time   = $(time.end_time)")
  println(io, "  Time step  = $(time.Δt)")
end

# exports
export ConstantTimeStepper

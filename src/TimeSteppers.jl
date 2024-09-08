"""
$(TYPEDEF)
"""
abstract type TimeStepper{T} end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct EndTimeException <: Exception
end

"""
$(TYPEDSIGNATURES)
"""
end_time_error() = throw(EndTimeException())

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct ConstantTimeStepper{T <: Number} <: TimeStepper{T}
  start_time::T
  end_time::T
  current_time::T
  current_time_step::Int
  Δt::T
end

"""
$(TYPEDSIGNATURES)
"""
function Base.similar(::ConstantTimeStepper)
  return ConstantTimeStepper(0.0, 0.0, 0.0, 0, 0.0)
end

"""
$(TYPEDSIGNATURES)
"""
function reset!(time::ConstantTimeStepper)
  time.current_time = time.start_time
  time.current_time_step = 1
end

"""
$(TYPEDSIGNATURES)
"""
function step!(time::ConstantTimeStepper) 
  temp = time.current_time + time.Δt
  time.current_time = temp
  time.current_time_step += 1
end

"""
$(TYPEDSIGNATURES)
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

abstract type TimeStepper{T} end

struct EndTimeException <: Exception
end

end_time_error() = throw(EndTimeException())

mutable struct ConstantTimeStepper{T <: Number} <: TimeStepper{T}
  start_time::T
  end_time::T
  current_time::T
  Δt::T
end

function reset!(time::ConstantTimeStepper)
  time.current_time = time.start_time
end

function step!(time::ConstantTimeStepper) 
  temp = time.current_time + time.Δt
  # if temp > time.end_time
  #   # error out here
  #   # end_time_error()
  #   @info "end time"
  # else
  time.current_time = temp
  # end
end

function ConstantTimeStepper(start_time::T, end_time::T, Δt::T) where T <: Number
  return ConstantTimeStepper(start_time, end_time, start_time, Δt)
end

function ConstantTimeStepper(input_settings::D) where D <: Dict
  @assert "start time" in keys(input_settings)
  @assert "end time"   in keys(input_settings)
  @assert "time step"  in keys(input_settings)
  
  start_time = input_settings["start time"]
  end_time   = input_settings["end time"]
  Δt         = input_settings["time step"]

  return ConstantTimeStepper(start_time, end_time, Δt)
end
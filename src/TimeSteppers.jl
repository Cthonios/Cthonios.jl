abstract type TimeStepper{T} end

struct EndTimeException <: Exception
end

end_time_error() = throw(EndTimeException())

mutable struct ConstantTimeStepper{T <: Number} <: TimeStepper{T}
  start_time::T
  end_time::T
  current_time::T
  current_time_step::Int
  Δt::T
end

function Base.similar(::ConstantTimeStepper)
  return ConstantTimeStepper(0.0, 0.0, 0.0, 0, 0.0)
end

function reset!(time::ConstantTimeStepper)
  time.current_time = time.start_time
  time.current_time_step = 1
end

function step!(time::ConstantTimeStepper) 
  temp = time.current_time + time.Δt
  # if temp > time.end_time
  #   # error out here
  #   # end_time_error()
  #   @info "end time"
  # else
  time.current_time = temp
  time.current_time_step += 1
  # end
end

function ConstantTimeStepper(start_time::T, end_time::T, Δt::T) where T <: Number
  return ConstantTimeStepper(start_time, end_time, start_time, 1, Δt)
end

function ConstantTimeStepper(input_settings::D) where D <: Dict  
  start_time = input_settings[Symbol("start time")]
  end_time   = input_settings[Symbol("end time")]
  Δt         = input_settings[Symbol("time step")]

  return ConstantTimeStepper(start_time, end_time, Δt)
end
abstract type AbstractTimeIntegrator end

# integrator methods
function Base.show(io::IO, int::T) where T <: AbstractTimeIntegrator
  println(io, "$T:")
  println(io, "  Start time = $(int.start_time)")
  println(io, "  End time   = $(int.end_time)")
  println(io, "  Time step  = $(int.Δt)")
end

current_time(int::AbstractTimeIntegrator) = int.current_time[1]
end_time(int::AbstractTimeIntegrator) = int.end_time
start_time(int::AbstractTimeIntegrator) = int.start_time
time_step(int::AbstractTimeIntegrator) = int.Δt[1]

function step!(integrator, solver, Uu, p)
  integration_step_header(integrator)
  step_new!(p, solver.objective)
  solve!(solver, Uu, p)
  return nothing
end

abstract type AbstractTimeIntegratorCache end

# cache methods
current_solution(cache::AbstractTimeIntegratorCache) = cache.U

include("Newmark.jl")
include("QuasiStatic.jl")

# exports
export Newmark
export QuasiStatic

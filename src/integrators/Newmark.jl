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

function step!(integrator::Newmark, solver, Uu, p)
  integration_step_header(integrator)

  # time updates
  temp = current_time(integrator) + time_step(integrator)
  integrator.current_time[1] = temp
  integrator.current_time_step[1] += 1

  # bc updates
  update_dirichlet_vals!(p, solver.objective)
  update_neumann_vals!(p, solver.objective)

  # predict
  predict!(p, integrator)

  # solve
  solve!(solver, Uu, p)

  # correct
  correct!(p, integrator)
  return nothing
end


# cache
struct NewmarkCache{Us, Uus}
  U::Us
  Uu_predicted::Uus
  U_predicted::Us
  Vu::Uus
  Au::Uus
end

function NewmarkCache(objective)
  U = create_fields(objective.domain)
  Uu_predicted = create_unknowns(objective.domain)
  U_predicted = create_fields(objective.domain)
  Vu = create_unknowns(objective.domain)
  Au = create_unknowns(objective.domain)
  return NewmarkCache(U, Uu_predicted, U_predicted, Vu, Au)
end

# newmark method specific methods
# def correct(UCorrection, V, A, dt):
#   A = UCorrection/(newmarkParameters.beta*dt*dt)
#   V += dt*newmarkParameters.gamma*A
#   return V, A

# TODO make non-allocating
# function correct(int::Newmark, Uu, Vu, Au)
function correct!(p, int::Newmark)
  @unpack β, γ, Δt = int
  # Vu_corr = copy(Vu)
  @unpack Vu, Au = p.integrator_cache
  Au .= Uu / (β * Δt * Δt)
  Vu += Δt * γ * Au
  return nothing
end

# def predict(U, V, A, dt):
  # U += dt*V + 0.5*dt*dt*(1.0 - 2.0*newmarkParameters.beta)*A
  # V += dt*(1.0 - newmarkParameters.gamma)*A
  # return U, V


# TODO make non-allocating
# function predict(int::Newmark, Uu, Vu, Au)
function predict!(p, int::Newmark)
  @unpack β, γ, Δt = int
  # Uu_pred = copy(Uu)
  # Vu_pred = copy(Vu)
  Uu_pred = p.integrator_cache.Uu_predicted
  Vu = p.integrator_cache.Vu

  Uu_pred += Δt * Vv + 0.5 * Δt * Δt * (1. - 2. * β) * Au
  Vu += Δt * (1. - γ) * Au
  return nothing
end

# connects types
integrator_cache(objective, ::Newmark) = NewmarkCache(objective)
integrator_unknowns(objective, ::Newmark) = (create_unknowns(objective.domain), create_unknowns(objective.domain))

Base.@kwdef struct ExplicitSolverSettings
    verbose::Bool = true
end

struct ExplicitSolver{O}
    objective_cache::O
    settings::ExplicitSolverSettings
end

function ExplicitSolver(objective_cache, p)
    return ExplicitSolver(objective_cache, ExplicitSolverSettings())
end

function solve!(solver::ExplicitSolver, u, p)
    Δt = p.times.Δt[1]
    o = solver.objective_cache
    v_old = o.solution_rate_old
    v = o.solution_rate
    a = o.solution_rate_rate
    M = o.lumped_mass

    internal_energy = _internal_energy(o, u, p)
    kinetic_energy = 0.5 * mapreduce(i -> M[i] * v[i]^2, +, eachindex(v))
    # TODO do external force (u dot f_ext is the energy)
    fill!(o.internal_energy, internal_energy)
    fill!(o.kinetic_energy, kinetic_energy)

    f_i = _internal_force(o, u, p)
    # TODO external forces
    a.data .= -f_i.data ./ M
    @. v = v_old + Δt * a
    @. u += Δt * v

    copyto!(v_old, v)
    return nothing
end


# function predict!(solver::ExplicitSolver, u, p)
#     # unpack a bunch of stuff
#     Δt = p.times.Δt[1]
#     γ = solver.objective_cache.γ
#     v = solver.objective_cache.solution_rate
#     a = solver.objective_cache.solution_rate_rate

#     u .= u .+= Δt * v .+ 0.5 * Δt * Δt * a
#     v .= v .+= (1. - γ) * Δt * a
    
#     # ensure bcs are correct
#     FiniteElementContainers.update_field_dirichlet_bcs!(u_pre, v_pre, p.dirichlet_bcs)
#     return nothing
# end

# function _update_acceleration(solver::ExplicitSolver, u, p)
#     # solver.objective_cache.inertial_force .= 
#     # unpack stuff
#     o = solver.objective_cache
#     v = o.solution_rate
#     a = o.solution_rate_rate
#     M = o.lumped_mass
    
#     internal_energy = _internal_energy(o, u, p)
#     kinetic_energy = 0.5 * dot(M, v .* v)
#     energy = internal_energy + kinetic_energy
#     # TODO do external force (u dot f_ext is the energy)
#     fill!(o.internal_energy, internal_energy)
#     fill!(o.kinetic_energy, kinetic_energy)

#     f_i = _internal_force(o, u, p)
#     f_k = M .* a
#     # TODO need external force
#     o.gradient .= f_i .+ f_k

    
# end

# function correct!(solver::ExplicitSolver, u, p)
#     # unpack a bunch of stuff
#     Δt = p.times.Δt[1]
#     γ = solver.objective_cache.γ
#     v = solver.objective_cache.solution_rate
#     a = solver.objective_cache.solution_rate_rate

#     # need to "solve" for acceleration

#     # a .= (u - u_pre) / β / Δt / Δt
#     v .+= γ * Δt * a

#     # ensure bcs are correct
#     FiniteElementContainers.update_field_dirichlet_bcs!(u, v, a, p.dirichlet_bcs)
# end

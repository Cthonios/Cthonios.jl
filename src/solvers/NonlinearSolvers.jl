include("linear_solvers/LinearSolvers.jl")

"""
"""
abstract type NonlinearSolverSettings end
"""
"""
abstract type NonlinearSolver{S, L} end

function logger end # To be defined for each solver
function solve! end # TO be defined for each solver
function update_unknown_dofs! end # To be defined for each solver

"""
$(TYPEDSIGNATURES)
"""
function update_unknown_dofs!(solver::NonlinearSolver, d::QuasiStaticDomain)
  # update the dofs using the linear solver method
  update_unknown_dofs!(solver.linear_solver, d)
end

# Wrappers methods around different domain types
"""
$(TYPEDSIGNATURES)
"""
function objective(solver, domain::QuasiStaticDomain, common, u)
  @timeit timer(common) "Objective" begin
    strain_energy!(solver, domain, domain.domain_cache, u, backend(common))
    o = domain.domain_cache.Π[1]
  end
  return o
end

"""
$(TYPEDSIGNATURES)
"""
function gradient(solver, domain::QuasiStaticDomain, common, u)
  @timeit timer(common) "Objective and gradient" begin
    internal_force!(solver, domain, domain.domain_cache, u, backend(common))
    g = @views domain.domain_cache.f[domain.dof.unknown_dofs]
  end
  return g
end

"""
$(TYPEDSIGNATURES)
"""
function hvp(solver, domain::QuasiStaticDomain, common, u, v)
  @timeit timer(common) "Hvp" begin
    stiffness_action!(solver, domain, domain.domain_cache, u, v, backend(common))
    Hv = @views domain.domain_cache.Hv[domain.dof.unknown_dofs]
  end
  return Hv
end

"""
$(TYPEDSIGNATURES)
"""
function hessian(solver, domain::QuasiStaticDomain, common, u)
  @timeit timer(common) "Hessian" begin
    stiffness!(solver, domain, domain.domain_cache, u, backend(common))
    K = SparseArrays.sparse!(solver.linear_solver.assembler) #|> Symmetric
  end
  return K
end

"""
$(TYPEDSIGNATURES)
"""
function objective_and_gradient(solver, domain::QuasiStaticDomain, common, u)
  @timeit timer(common) "Objective and gradient" begin
    strain_energy_and_internal_force!(solver, domain, domain.domain_cache, u, backend(common))
    o = domain.domain_cache.Π[1]
    g = @views domain.domain_cache.f[domain.dof.unknown_dofs]
  end
  return o, g
end

"""
$(TYPEDSIGNATURES)
"""
function objective_gradient_and_hessian(solver, domain::QuasiStaticDomain, common, u)
  @timeit timer(common) "Objective, gradient, and Hessian" begin
    strain_energy_internal_force_and_stiffness!(solver, domain, domain.domain_cache, u, backend(common))
    o = domain.domain_cache.Π[1]
    g = @views domain.domain_cache.f[domain.dof.unknown_dofs]
    K = SparseArrays.sparse!(solver.linear_solver.assembler)
  end
  return o, g, K
end

include("NewtonSolver.jl")
include("TrustRegionSolver.jl")

function setup_nonlinear_solver(
  inputs::D, domain::QuasiStaticDomain, backend
) where {D <: Dict{Symbol, Any}}
  type = eval(Meta.parse(inputs[:type]))
  return type(inputs, domain, backend)
end

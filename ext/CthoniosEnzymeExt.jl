module CthoniosEnzymeExt

using Cthonios
using Enzyme
using FiniteElementContainers
using LinearAlgebra
using LinearSolve
using UnPack

# was using Enzyme 11.20

# helpers
function setup_ddomain(::ReverseMode, domain::Cthonios.Domain)
  ddomain = similar(domain)
  ddomain.X .= zero(eltype(ddomain.X))
  ddomain.U .= zero(eltype(ddomain.U))
  ddomain.U_bc .= zero(eltype(ddomain.U_bc))
  ddomain.props .= zero(eltype(ddomain.props))
  # ddomain.state_old .= zero(eltype(ddomain.state_old))
  # ddomain.state_new .= zero(eltype(ddomain.state_new))
  ddomain.solver_cache.Π .= one(eltype(ddomain.solver_cache.Π))
  ddomain.solver_cache.Πs .= zero(eltype(ddomain.solver_cache.Πs))
  ddomain.solver_cache.assembler.residuals .= 
    zero(eltype(ddomain.solver_cache.assembler.residuals))
  ddomain.solver_cache.assembler.stiffnesses .= 
    zero(eltype(ddomain.solver_cache.assembler.stiffnesses))
  return ddomain
end

function setup_dUu(::ForwardMode, Uu)
  dUu = similar(Uu)
  dUu .= one(eltype(dUu))
  return dUu
end

function setup_dUu(::ReverseMode, Uu)
  dUu = similar(Uu)
  dUu .= zero(eltype(dUu))
  return dUu
end

# first order methods
function Cthonios.gradient!(domain, ddomain, Uu, dUu)
  autodiff_deferred(
    Reverse, Cthonios.internal_energy!,
    Duplicated(domain, ddomain),
    Duplicated(Uu, dUu)
  )
  return nothing
end

function Cthonios.gradient(mode::ReverseMode, domain::Domain, Uu)
  ddomain = setup_ddomain(mode, domain)
  dUu = setup_dUu(mode, Uu)
  Cthonios.gradient!(domain, ddomain, Uu, dUu)
  return ddomain, dUu
end

function Cthonios.gradient(mode::ReverseMode, prob::P, Uu) where P <: Cthonios.AbstractProblem
  Cthonios.gradient(mode, prob.domain, Uu)
end

# attempt at AD over single nonlinear step
function grad_over_solve!(Uu, dUu, objective, solver, dsolver, domain, ddomain)
  autodiff(
    Reverse, Cthonios.solve!,
    Duplicated(solver, dsolver),
    Const(objective),
    Duplicated(domain, ddomain),
    Duplicated(Uu, dUu)
    # linear_solver_alg = Const(CHOLMODFactorization())
  )
  return nothing
end 

function Cthonios.solve(mode::ReverseMode, prob::P, Uu) where P <: Cthonios.AbstractProblem
  @unpack solver, objective, domain = prob
  solver = solver.linear_solver
  Cthonios.step!(domain)
  ddomain = setup_ddomain(mode, domain)
  dUu = setup_dUu(mode, Uu)
  dsolver = similar(solver)
  grad_over_solve!(Uu, dUu, objective, solver, dsolver, domain, ddomain)
  return ddomain, dUu
end

# TOOD LLVM barf for all below
# second order methods
function residual_dot_v!(r_dot_v, domain, ddomain, Uu, dUu, Vv)
  Cthonios.gradient!(domain, ddomain, Uu, dUu)
  r_dot_v[1] = dot(domain.solver_cache.assembler.residuals[domain.dof.unknown_dofs], Vv)
  return nothing
end

function Cthonios.hvp(::ReverseMode, prob::P, Uu, Vv) where P <: Cthonios.AbstractProblem
  # @unpack static, cache = prob.domain
  domain = prob.domain
  # dcache = setup_dcache(Reverse, cache)
  bdomain = setup_ddomain(Reverse, domain)
  bUu = setup_dUu(Reverse, Uu)
  dUu = setup_dUu(Forward, Uu)
  dbUu = setup_dUu(Reverse, Uu)
  r_dot_v = zeros(Float64, 1)
  dr_dot_v = ones(Float64, 1)
  autodiff(
    # Forward, residual_dot_v!,
    # Duplicated(r_dot_v, dr_dot_v),
    Forward, Cthonios.gradient!,
    Const(domain),
    Const(bdomain),
    Duplicated(Uu, dUu),
    Duplicated(bUu, Vv),
    # Const(Vv)
  )
  # dcache, dUu
end

# function Cthonios.residual!(dcache, dUu, ::ReverseMode, static, cache, Uu)
#   autodiff(
#     Reverse, Cthonios.internal_energy!, 
#     Const(static), 
#     Duplicated(cache, dcache),
#     Duplicated(Uu, dUu)
#   )
#   return nothing
# end

# function Cthonios.residual(::ReverseMode, prob::P, Uu) where P <: Cthonios.AbstractProblem
#   @unpack static, cache = prob.domain
#   dcache = setup_dcache(Reverse, cache)
#   dUu = setup_dUu(Reverse, Uu)
#   Cthonios.residual!(dcache, dUu, Reverse, static, cache, Uu)
#   return dUu
# end

# # general methods helpers
# # function grad!(dcache, dUu, Reverse, obj, static, cache, Uu)
# function grad!(obj, static, cache, dcache, Uu, dUu)
#   autodiff_deferred(
#     Reverse, obj.value,
#     Const(static),
#     Duplicated(cache, dcache),
#     Duplicated(Uu, dUu)
#   )
#   return nothing
# end

# # function residual_dot_v!(r_dot_v, static, cache, Uu, Vv)
#   # Cthonios.internal_force!(static, cache, Uu)
# # function grad_dot_v!(r_dot_v, dcache, dUu, ::ReverseMode, obj, static, cache, Uu, Vv)
# function grad_dot_v!(obj, static, r_dot_v, cache, dcache, Uu, dUu, Vv)
#   # grad!(dcache, dUu, Reverse, obj, static, cache, Uu)
#   grad!(obj, static, cache, dcache, Uu, dUu)
#   # r_dot_v[1] = dot(cache.solver_cache.assembler.residuals[static.dof.unknown_dofs], Vv)
#   r_dot_v[1] = dot(dUu, Vv)
#   return nothing
# end

# function Cthonios.grad(::ReverseMode, prob::P, Uu) where P <: Cthonios.AbstractProblem
#   @unpack static, cache = prob.domain
#   dcache = setup_dcache(Reverse, cache)
#   dUu = setup_dUu(Reverse, Uu)
#   grad!(prob.objective, static, cache, dcache, Uu, dUu)
#   dcache, dUu
# end

# # fix below to not have reverse
# function Cthonios.hvp(::ReverseMode, prob::P, Uu, Vv) where P <: Cthonios.AbstractProblem
#   @unpack static, cache = prob.domain

#   bcache = setup_dcache(Reverse, cache)
#   dcache = setup_dcache(Reverse, cache)
#   dcache.solver_cache.Π .= zeros(Float64, 1)
#   dbcache = setup_dcache(Reverse, cache)
#   dbcache.solver_cache.Π .= zeros(Float64, 1)

#   bUu = setup_dUu(Reverse, Uu)
#   dUu = setup_dUu(Forward, Uu)
#   dbUu = setup_dUu(Reverse, Uu)

#   r_dot_v = zeros(Float64, 1)
#   dr_dot_v = ones(Float64, 1)

#   # puts the grad in bcache and bUu
#   grad_dot_v!(prob.objective, static, r_dot_v, cache, bcache, Uu, bUu, Vv)
  
#   autodiff(
#     Forward, grad_dot_v!,
#     Const(prob.objective),
#     Const(static),
#     Duplicated(r_dot_v, dr_dot_v),
#     # Duplicated(cache, dcache),
#     # Duplicated(bcache, dbcache),
#     Const(cache),
#     Const(bcache),
#     Duplicated(Uu, dUu),
#     Duplicated(bUu, dbUu),
#     Const(Vv)
#   )
#   # autodiff(
#   #   Forward, grad_dot_v!,
#   #   Const(static),
#   #   Duplicated(r_dot_v, dr_dot_v),
#   #   Duplicated(x, )
#   #   # Duplicated()
#   #   # Duplicated(Duplicated(cache, bcache), Duplicated(dcache, dbcache)),
#   #   # Duplicated(Duplicated(Uu, bUu), Duplicated(dUu, dbUu))
#   # )
# end

end # module
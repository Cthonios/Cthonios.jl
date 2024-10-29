module CthoniosEnzymeExt

using Cthonios
using DifferentiationInterface
using Enzyme
using LinearAlgebra
using Parameters
using SparseArrays


function gradient_deferred!(
  dval, dUu, dp,
  val, o, Uu, p
)
  autodiff_deferred(
    Reverse, Cthonios.objective!,
    Duplicated(val, dval),
    Const(o),
    Duplicated(Uu, dUu),
    Duplicated(p, dp)
  )
  return nothing
end

function Cthonios.grad_u(o::Objective, Uu, p)
  val, dval = zeros(1), zeros(1)
  dUu, dp = make_zero(Uu), make_zero(p)
  gradient_deferred!(dval, dUu, dp, val, o, Uu, p)
  dUu
end

function Cthonios.grad_p(o::Objective, Uu, p)
  val, dval = zeros(1), zeros(1)
  dUu, dp = make_zero(Uu), make_zero(p)
  gradient_deferred!(dval, dUu, dp, val, o, Uu, p)
  dp
end

function Cthonios.hvp_u(o::Objective, Uu, p, Vv)
  val, dval = zeros(1), zeros(1)
  grad_u, grad_p = make_zero(Uu), make_zero(p)
  autodiff(
    Forward, gradient_deferred!, Const(Reverse),
    Duplicated(),
    Const(o),
    Duplicated(),
    Duplicated()
  )
end

# out of place objective methods
# function Cthonios.grad_u(o::Objective, Uu, p)
# # function Cthonios.grad_u(f, domain, Uu, p)
#   backend = AutoEnzyme(Reverse)
#   func = x -> Cthonios.objective(o, x, p)
#   # func = x -> Cthonios.domain_iterator(f, domain, x, p)
#   DifferentiationInterface.gradient(func, backend, Uu)
# end


# function gradient_deferred!(::ReverseMode, dΠ, dUu, dp, objective, Π, Uu, p)
#   autodiff_deferred(
#     Reverse, Cthonios.objective!, 
#     Duplicated(Π, dΠ),
#     Const(objective),
#     Duplicated(Uu, dUu),
#     Duplicated(p, dp)
#   )
#   return nothing
# end

# function Cthonios.gradient(::ReverseMode, solver, Uu, p)
#   dUu = make_zero(Uu)
#   dp = make_zero(p)
#   Π = solver.o
#   dΠ = make_zero(Π)
#   dΠ .= one(eltype(dΠ))
#   gradient_deferred!(Reverse, dΠ, dUu, dp, solver.objective, Π, Uu, p)
#   return dUu, dp
# end

function Cthonios.solve!(warm_start::Cthonios.WarmStart, solver, objective, Uu, p)
  @info "Warm start"
  @unpack Uu_old, U_old = warm_start
  Uu_old .= zero(eltype(Uu_old))
  U_old .= zero(eltype(U_old))
  U = p.U

  # field unknowns, same for both
  Cthonios.update_field_unknowns!(U_old, objective.domain, Uu)
  Cthonios.update_field_unknowns!(U, objective.domain, Uu)

  Cthonios.update_field_bcs!(U_old, objective.domain, p.Ubc)
  Cthonios.step!(p)
  Cthonios.update_dirichlet_vals!(p, objective)
  Cthonios.update_field_bcs!(U, objective.domain, p.Ubc)

  # get bc increment
  dU = similar(U)
  dU .= zero(eltype(dU))
  @. dU = U_old - p.U

  # need stiffness array
  R = solver.assembler.residuals
  Cthonios.update_preconditioner!(solver, objective, Uu, p)

  # set up seed arrays
  dR = make_zero(R)
  dp = make_zero(p)

  # run enzyme
  autodiff(
    Reverse, Cthonios.gradient!,
    Duplicated(R, dR), 
    Const(objective),
    # Duplicated(U, dU),
    Duplicated(p, dp)
    # Duplicated(U, dU),
    # Const(p.X)
  )

  # a little ugly below
  solver.assembler.stiffnesses .= zero(eltype(solver.assembler.stiffnesses))
  Cthonios.domain_iterator!(solver.assembler, objective.hessian, objective.domain, Uu, p)
  K = SparseArrays.sparse!(solver.assembler) |> Symmetric

  b = R[objective.domain.dof.unknown_dofs]
  temp_Uu = K \ b

  @. Uu -= temp_Uu
  return nothing
end

end # module

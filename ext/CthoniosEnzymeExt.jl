module CthoniosEnzymeExt

using Cthonios
using Enzyme
using LinearAlgebra
using Parameters
using SparseArrays

function gradient_deferred!(::ReverseMode, dΠ, dUu, dp, objective, Π, Uu, p)
  autodiff_deferred(
    Reverse, Cthonios.objective!, 
    Duplicated(Π, dΠ),
    Const(objective),
    Duplicated(Uu, dUu),
    Duplicated(p, dp)
  )
  return nothing
end

function Cthonios.gradient(::ReverseMode, solver, Uu, p)
  dUu = make_zero(Uu)
  dp = make_zero(p)
  Π = solver.o
  dΠ = make_zero(Π)
  dΠ .= one(eltype(dΠ))
  gradient_deferred!(Reverse, dΠ, dUu, dp, solver.objective, Π, Uu, p)
  return dUu, dp
end

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

  # run enzyme
  autodiff(
    Reverse, Cthonios.gradient!,
    Duplicated(R, dR), 
    Const(objective),
    Duplicated(U, dU),
    Const(p.X)
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

abstract type AbstractPreconditioner end

mutable struct CholeskyPreconditioner{A, P} <: AbstractPreconditioner
  assembler::A
  preconditioner::P
end

function CholeskyPreconditioner(obj::Objective)
  asm = StaticAssembler(obj.domain)
  # update_unknown_dofs!(obj.domain)
  update_unknown_dofs!(obj.domain, asm)
  Uu = create_unknowns(obj.domain)
  U = create_fields(obj.domain)
  # TODO need stuff here
  p = nothing
  update_fields!(U, obj.domain, Uu)
  domain_iterator!(asm, hessian, obj.domain, U, p)
  H = SparseArrays.sparse!(asm) |> Symmetric
  P = cholesky(H)
  return CholeskyPreconditioner(asm, P)
end

function LinearAlgebra.ldiv!(y, P::CholeskyPreconditioner, v)
  y .= P.preconditioner \ v
  return nothing
end

function update_preconditioner!(P::CholeskyPreconditioner, obj, Uu, p)
  P.assembler.stiffnesses .= zero(eltype(P.assembler.stiffnesses))
  U = create_fields(obj.domain)
  update_bcs!(U, obj.domain)
  update_fields!(U, obj.domain, Uu)
  domain_iterator!(P.assembler, hessian, obj.domain, U, p)
  H = SparseArrays.sparse!(P.assembler) |> Symmetric

  attempt = 1
  while attempt < 10
    @info "Updating preconditioner, attempt = $attempt"
    try
      if attempt == 1
        P.preconditioner = cholesky(H)
      else
        shift = 10.0^(-5 + attempt)
        P.preconditioner = cholesky(H; shift=shift)
      end
      return nothing
    catch e
      @info e
      @info "Failed to factor preconditioner. Attempting again"
      attempt += 1 
    end
  end
  return nothing
end


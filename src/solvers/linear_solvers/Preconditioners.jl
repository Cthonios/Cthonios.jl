abstract type AbstractPreconditioner end

mutable struct CholeskyPreconditioner{A, P} <: AbstractPreconditioner
  assembler::A
  preconditioner::P
end

function CholeskyPreconditioner(obj::Objective, p)
  asm = StaticAssembler(obj.domain)
  update_unknown_dofs!(obj.domain, asm)
  # TODO need stuff here
  # inefficiency here by creating these copies
  Uu = create_unknowns(obj.domain)
  H = hessian!(asm, obj, Uu, p)
  P = cholesky(H)
  return CholeskyPreconditioner(asm, P)
end

function LinearAlgebra.ldiv!(y, P::CholeskyPreconditioner, v)
  y .= P.preconditioner \ v
  return nothing
end

function update_preconditioner!(P::CholeskyPreconditioner, obj, Uu, p)
  H = hessian!(P.assembler, obj, Uu, p)
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

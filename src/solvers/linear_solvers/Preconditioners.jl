abstract type AbstractPreconditioner end

mutable struct CholeskyPreconditioner{A, P, T} <: AbstractPreconditioner
  assembler::A
  preconditioner::P
  timer::T
end

function CholeskyPreconditioner(obj::Objective, p, timer)
  @timeit timer "CholeskyPreconditioner - setup" begin
    asm = StaticAssembler(obj.domain)
    update_unknown_dofs!(obj.domain, asm)
    # TODO need stuff here
    # inefficiency here by creating these copies
    Uu = create_unknowns(obj.domain)
    H = hessian!(asm, obj, Uu, p)
    P = cholesky(H)
  end
  return CholeskyPreconditioner(asm, P, timer)
end

timer(P::CholeskyPreconditioner) = P.timer

function LinearAlgebra.ldiv!(y, P::CholeskyPreconditioner, v)
  @timeit timer(P) "CholeskyPreconditioner - ldiv!" begin
    y .= P.preconditioner \ v
  end
  return nothing
end

function update_preconditioner!(P::CholeskyPreconditioner, obj, Uu, p)
  @timeit timer(P) "CholeskyPreconditioner - update_preconditioner!" begin
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
  end
  return nothing
end

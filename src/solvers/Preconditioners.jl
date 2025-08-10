abstract type AbstractPreconditioner end
function update_preconditioner! end

# fall back methods
LinearAlgebra.ldiv!(y, P::AbstractPreconditioner, x) = ldiv!(y, P.preconditioner, x)
timer(P::T) where T <: AbstractPreconditioner = P.timer

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct CholeskyPreconditioner{A, P, T} <: AbstractPreconditioner
  assembler::A
  preconditioner::P
  timer::T
end

"""
$(TYPEDSIGNATURES)
"""
function CholeskyPreconditioner(obj::AbstractObjectiveCache, p, timer)
  @timeit timer "CholeskyPreconditioner - setup" begin
    # asm = StaticAssembler(obj.domain)
    # asm = obj.assembler
    asm = assembler(obj)
    # update_unknown_dofs!(obj.domain, asm)
    # TODO need stuff here
    # inefficiency here by creating these copies
    Uu = create_unknowns(asm, H1Field)
    H = hessian(obj, Uu, p)
    P = cholesky(H)
  end
  return CholeskyPreconditioner(asm, P, timer)
end

"""
$(TYPEDSIGNATURES)
"""
function LinearAlgebra.ldiv!(y, P::CholeskyPreconditioner, v)
  @timeit timer(P) "CholeskyPreconditioner - ldiv!" begin
    y .= P.preconditioner \ v
  end
  return nothing
end

function update_preconditioner!(P::CholeskyPreconditioner, obj, Uu, p; verbose=false)
  @timeit timer(P) "CholeskyPreconditioner - update_preconditioner!" begin
    # H = hessian!(P.assembler, obj, Uu, p)
    H = hessian(obj, Uu, p)
    attempt = 1
    while attempt < 10
      if verbose
        @info "Updating preconditioner, attempt = $attempt"
      end
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

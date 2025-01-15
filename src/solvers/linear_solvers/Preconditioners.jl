"""
$(TYPEDEF)
Abstract base preconditioner type.
"""
abstract type AbstractPreconditioner end

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
function CholeskyPreconditioner(obj::AbstractObjective, p, timer)
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

"""
$(TYPEDSIGNATURES)
"""
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

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct LimitedLDLPreconditioner{A, P, T} <: AbstractPreconditioner
  assembler::A
  preconditioner::P
  timer::T
end

"""
$(TYPEDSIGNATURES)
"""
function LimitedLDLPreconditioner(obj::AbstractObjective, p, timer)
  @timeit timer "CholeskyPreconditioner - setup" begin
    asm = StaticAssembler(obj.domain)
    update_unknown_dofs!(obj.domain, asm)
    # TODO need stuff here
    # inefficiency here by creating these copies
    Uu = create_unknowns(obj.domain)
    hessian!(asm, obj, Uu, p)
    H = SparseArrays.sparse!(asm)
    P = lldl(H; α=1e-5, α_increase_factor=10)
  end
  return LimitedLDLPreconditioner(asm, P, timer)
end

function update_preconditioner!(P::LimitedLDLPreconditioner, obj, Uu, p)
  @timeit timer(P) "CholeskyPreconditioner - update_preconditioner!" begin
    hessian!(P.assembler, obj, Uu, p)
    H = SparseArrays.sparse!(P.assembler) #|> Symmetric
    attempt = 1
    @info "Updating preconditioner, attempt = $attempt"
    P.preconditioner = lldl(H; α=1e-5, α_increase_factor=10)
  end
  return nothing
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct LDLPreconditioner{A, P, T} <: AbstractPreconditioner
  assembler::A
  preconditioner::P
  timer::T
end

"""
$(TYPEDSIGNATURES)
"""
function LDLPreconditioner(obj::AbstractObjective, p, timer)
  @timeit timer "CholeskyPreconditioner - setup" begin
    asm = StaticAssembler(obj.domain)
    update_unknown_dofs!(obj.domain, asm)
    # TODO need stuff here
    # inefficiency here by creating these copies
    Uu = create_unknowns(obj.domain)
    H = hessian!(asm, obj, Uu, p)
    P = ldl(SparseArrays.sparse!(asm))
  end
  return LDLPreconditioner(asm, P, timer)
end

# TODO figure out how to apply shift here
function update_preconditioner!(P::LDLPreconditioner, obj, Uu, p)
  @timeit timer(P) "LDLPreconditioner - update_preconditioner!" begin
    H = hessian!(P.assembler, obj, Uu, p)
    ldl_factorize!(SparseArrays.sparse!(P.assembler), P.preconditioner)
  end
  return nothing
end

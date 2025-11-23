abstract type AbstractPreconditioner end
function update_preconditioner! end

# fall back methods
LinearAlgebra.ldiv!(y, P::AbstractPreconditioner, x) = ldiv!(y, P.preconditioner, x)
timer(P::T) where T <: AbstractPreconditioner = P.timer

struct NoPreconditioner <: AbstractPreconditioner
  timer::TimerOutput
end

function NoPreconditioner(obj, timer)
  return NoPreconditioner(timer)
end

function LinearAlgebra.ldiv!(y, P::NoPreconditioner, v)
  @timeit timer(P) "NoPreconditioner - ldiv!" begin
    ldiv!(y, I, v)
  end
  return nothing
end

function update_preconditioner!(P::NoPreconditioner, obj, Uu, p; verbose=false)
  return nothing
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct CholeskyPreconditioner{A, P, T} <: AbstractPreconditioner
  assembler::A
  preconditioner::P
  timer::T
end

function _cholesky(A::SparseMatrixCSC; shift=0.0)
  return cholesky(A |> Symmetric; shift=shift)
end

function _cholesky!(F, A::SparseMatrixCSC; shift=0.0)
  cholesky!(F, A |> Symmetric; shift=shift)
end

"""
$(TYPEDSIGNATURES)
"""
function CholeskyPreconditioner(obj, timer)
  @timeit timer "CholeskyPreconditioner - setup" begin
    asm = assembler(obj)
    p = parameters(obj)
    # update_unknown_dofs!(obj.domain, asm)
    # TODO need stuff here
    # inefficiency here by creating these copies
    Uu = create_unknowns(asm)
    # H = hessian(obj, Uu, p)
    # P = cholesky(H)
    assemble_stiffness!(asm, obj.objective.hessian_u, Uu, p)
    H = stiffness(asm)
    P = _cholesky(H)
  end
  return CholeskyPreconditioner(asm, P, timer)
end

"""
$(TYPEDSIGNATURES)
"""
function LinearAlgebra.ldiv!(y, P::CholeskyPreconditioner, v)
  @timeit timer(P) "CholeskyPreconditioner - ldiv!" begin
    # y .= P.preconditioner \ v
    _ldiv!(y, P, v, KA.get_backend(y))
  end
  return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function _ldiv!(y, P::CholeskyPreconditioner, v, ::CPU)
  y .= P.preconditioner \ v

  # below doesn't work. I guess we just have to accept these allocations
  # ldiv!(y, P.preconditioner.L, v)
  # ldiv!(y, P.preconditioner.L')
  return nothing
end

function update_preconditioner!(P::CholeskyPreconditioner, obj, Uu, p; verbose=false)
  @timeit timer(P) "CholeskyPreconditioner - update_preconditioner!" begin
    # H = hessian!(P.assembler, obj, Uu, p)
    # H = hessian(obj, Uu, p)
    asm = assembler(obj)
    assemble_stiffness!(asm, obj.objective.hessian_u, Uu, p)
    H = stiffness(asm)
    attempt = 1
    while attempt < 10
      if verbose
        @info "Updating preconditioner, attempt = $attempt"
      end
      try
        if attempt == 1
          # P.preconditioner = _cholesky(H)
          _cholesky!(P.preconditioner, H)
        else
          shift = 10.0^(-5 + attempt)
          # P.preconditioner = _cholesky(H; shift=shift)
          _cholesky!(P.preconditioner, H; shift=shift)
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

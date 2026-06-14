abstract type AbstractPreconditioner end
function update_preconditioner! end

# fall back methods
LinearAlgebra.ldiv!(y, P::AbstractPreconditioner, x) = ldiv!(y, P.preconditioner, x)
timer(P::T) where T <: AbstractPreconditioner = P.timer

# mutable struct AMGPreconditioner{P} <: AbstractPreconditioner
#   preconditioner::P
#   timer::TimerOutput
# end

# function AMGPreconditioner(obj::AbstractObjective, u, p, timer = TimerOutput())
#   H = hessian(obj, u, p)
#   ml = ruge_stuben(H)
#   return AMGPreconditioner(aspreconditioner(ml), timer)
# end

# function LinearAlgebra.ldiv!(y, P::AMGPreconditioner, x)
#   ldiv!(y, P.preconditioner, x)
#   return nothing
# end

# function update_preconditioner!(P::AMGPreconditioner, obj, u, p; verbose = false)
#   H = hessian(obj, u, p)
#   ml = ruge_stuben(H)
#   P.preconditioner = aspreconditioner(ml)
#   return nothing
# end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct CholeskyPreconditioner{P} <: AbstractPreconditioner
  preconditioner::P
  timer::TimerOutput
end

"""
$(TYPEDSIGNATURES)
"""
function CholeskyPreconditioner(obj::AbstractObjective, u, p, timer = TimerOutput())
  @timeit timer "CholeskyPreconditioner - setup" begin
    H = hessian(obj, u, p) |> Symmetric
    P = cholesky(H)
  end
  return CholeskyPreconditioner{typeof(P)}(P, timer)
end

"""
$(TYPEDSIGNATURES)
"""
function LinearAlgebra.ldiv!(y, P::CholeskyPreconditioner, v)
  y .= P.preconditioner \ v

  # below doesn't work. I guess we just have to accept these allocations
  # ldiv!(y, P.preconditioner.L, v)
  # ldiv!(y, P.preconditioner.L')
  return nothing
end

function update_preconditioner!(P::CholeskyPreconditioner, obj, u, p; verbose = false)
  @timeit timer(P) "CholeskyPreconditioner - update_preconditioner!" begin
    H = hessian(obj, u, p) |> Symmetric
    attempt = 1
    while attempt < 10
      if verbose
        @info "Updating preconditioner, attempt = $attempt"
      end
      try
        if attempt == 1
          cholesky!(P.preconditioner, H)
        else
          shift = 10.0^(-5 + attempt)
          # cholesky!(P.preconditioner, H; shift=shift)
          cholesky!(P.preconditioner, H + Diagonal(shift .* abs(diag(H))))
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

struct JacobiPreconditioner{P} <: AbstractPreconditioner
  inv_diag::P
  timer::TimerOutput
end

# need inplace diagonal method
function JacobiPreconditioner(objective, u, p, timer = TimerOutput())
  FEC.assemble_diagonal!(assembler(objective), mass!, u, p)
  d_mass = FEC.diagonal(asm)
  inv_diag = 1. ./ d_mass
  return JacobiPreconditioner(inv_diag, timer)
end

function LinearAlgebra.ldiv!(y, P::JacobiPreconditioner, x)
  @. y = P.inv_diag * x
  return nothing
end

function update_preconditioner!(P::JacobiPreconditioner, objective, u, p; verbose = false)
  FEC.assemble_diagonal!(assembler(objective), mass, u, p)
  d = FEC.diagonal(asm)
  P.inv_diag .= 1. ./ d
  return nothing
end

struct LLDLPreconditioner{P} <: AbstractPreconditioner
  preconditioner::P
  timer::TimerOutput
end

function LLDLPreconditioner(objective, u, p, timer = TimerOutput())
  H = hessian_symmetric(objective, u, p)
  P = LimitedLDLFactorizations.lldl(H)
  return LLDLPreconditioner(P, timer)
end

function update_preconditioner!(P::LLDLPreconditioner, objective, u, p; verbose = false)
  @timeit timer(P) "LLDLPreconditioner - update_preconditioner!" begin
    H = hessian_symmetric(objective, u, p)
    attempt = 1
    while attempt < 10
      if verbose
        @info "Updating preconditioner, attempt = $attempt"
      end
      try
        if attempt == 1
          lldl_factorize!(P.preconditioner, H)
        else
          shift = 10.0^(-5 + attempt)
          update_shift!(P.preconditioner, shift)
          lldl_factorize!(P.preconditioner, H)
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
struct LUPreconditioner{P} <: AbstractPreconditioner
  preconditioner::P
  timer::TimerOutput
end

"""
$(TYPEDSIGNATURES)
"""
function LUPreconditioner(obj::AbstractObjective, u, p, timer = TimerOutput())
  @timeit timer "LUPreconditioner - setup" begin
    H = hessian(obj, u, p)
    P = lu(H)
  end
  return LUPreconditioner{typeof(P)}(P, timer)
end

"""
$(TYPEDSIGNATURES)
"""
function LinearAlgebra.ldiv!(y, P::LUPreconditioner, v)
  y .= P.preconditioner \ v
  return nothing
end


function update_preconditioner!(P::LUPreconditioner, objective, u, p; verbose = false)
  @timeit timer(P) "LUPreconditioner - update_preconditioner!" begin
    H = hessian(objective, u, p)
    lu!(P.preconditioner, H)
  end
  return nothing
end

struct NoPreconditioner <: AbstractPreconditioner
  timer::TimerOutput
end

function NoPreconditioner(obj, u, p, timer = TimerOutput())
  return NoPreconditioner(timer)
end

function LinearAlgebra.ldiv!(y, P::NoPreconditioner, v)
  ldiv!(y, I, v)
  return nothing
end

function update_preconditioner!(P::NoPreconditioner, obj, u, p; verbose=false)
  return nothing
end

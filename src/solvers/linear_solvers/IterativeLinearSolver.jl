# Preconditioner below is copy-pasta from LinearSolve
# LinearSolve depends upon LoopVectorizations which makes compiling
# to an executable difficult
# source: https://github.com/SciML/LinearSolve.jl/blob/main/src/preconditioners.jl
#
struct InvPreconditioner{T}
  P::T
end

Base.eltype(A::InvPreconditioner) = Base.eltype(A.P)
LinearAlgebra.ldiv!(A::InvPreconditioner, x) = mul!(x, A.P, x)
LinearAlgebra.ldiv!(y, A::InvPreconditioner, x) = mul!(y, A.P, x)
LinearAlgebra.mul!(y, A::InvPreconditioner, x) = ldiv!(y, A.P, x)

####################################################################

# My attempt at a linear operator
# this is highly specific to Cthonios

# struct DomainLinearOperator{Domain}

# end

# function LinearAlgebra.mul!(y, domain::QuasiStaticDomain, x)

# end

##################################################################

struct IterativeLinearSolverSettings <: AbstractLinearSolverSettings
end

struct IterativeLinearSolver{Settings, Assembler} <: AbstractLinearSolver{Settings, Assembler}
  settings::Settings
  assembler::Assembler
end

# TODO maybe add some settings
function IterativeLinearSolver(input_settings::Dict{Symbol, Any}, domain::QuasiStaticDomain)
  settings = IterativeLinearSolverSettings()
  assembler = MatrixFreeStaticAssembler(domain.dof)
  update_unknown_dofs!(domain) # need to update prior to making an operator

  # also need to assemble once
  # U = domain.domain_cache.U
  # V = FiniteElementContainers.create_fields(domain.dof)
  Uu = FiniteElementContainers.create_unknowns(domain.dof) # TODO remove this allocation
  Vv = FiniteElementContainers.create_unknowns(domain.dof)
  # state = domain.domain_cache.state
  # props = domain.domain_cache.props
  @unpack U, state, props, Π, V = domain.domain_cache

  # @assert false
  # stiffness_action!(assembler, U, domain, Uu, state, props, domain.coords, V, Vv)
  NU = length(domain.dof.unknown_dofs)
  op = LinearOperator(
    Float64, NU, NU, true, true,
    (r, v) -> stiffness_action!(solver, domain, Uu, Vv)
  )

  P = InvPreconditioner(op)
  Yy = FiniteElementContainers.create_unknowns(domain.dof)
  display(P)
  ldiv!(Yy, P, Vv)
  display(Yy)
  @assert false
end

mutable struct LinearSolver{
  Preconditioner,
  Residual,
  Stiffness,
  Unknowns
}
  precond::Preconditioner
  residual::Residual
  stiffness::Stiffness
  unknowns::Unknowns
end 

function LinearSolver(
  input_settings::D, domain::QuasiStaticDomain
) where D <: Dict

  Uu = create_unknowns(domain)
  R  = create_unknowns(domain)
  # TODO need to set up custom solvers
  # R = residual(domain, Uu)
  # K = stiffness(domain, Uu)
  K = stiffness(domain)

  # setup solver
  # if :type in keys(input_settings)
  #   if input_settings[:type] == "default"
  #     solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
  #   else
  #     solver = eval(Meta.parse(input_settings[:type]))()
  #   end
  # else
  #   solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
  # end

  # set up preconditioner
  if :preconditioner in keys(input_settings)
    if input_settings[:preconditioner] == "default"
      # Pl = IterativeSolvers.IdentityOperator(length(Uu))
      Pl = I
    elseif input_settings[:preconditioner] == "cholesky"
      Pl = cholesky(K)
    else
      Pl = eval(Meta.parse(input_settings[:preconditioner]))(K)
    end
  else
    # Pl = IdentityOperator(length(Uu))
    Pl = I
  end

  return LinearSolver(Pl, R, K, Uu)
end

function Base.show(io::IO, solver::LinearSolver)
  print(io, "LinearSolver\n", 
            "  Preconditioner = $(solver.precond)")
end

function update_residual!(
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}

  R = domain.assembler.residuals
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  internal_force!(R, U, domain, Uu, state, props, domain.coords)
  solver.residual .= @views R[domain.dof.unknown_dofs]
end

function update_stiffness!(
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}

  # stiffness!(domain, Uu, domain.coords, U)
  state = domain.domain_cache.state
  props = domain.domain_cache.props
  stiffness!(domain.assembler, U, domain, Uu, state, props, domain.coords)
  solver.stiffness = SparseArrays.sparse!(domain.assembler)
end

function update_residual_and_stiffness!(
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}

  state = domain.domain_cache.state
  props = domain.domain_cache.props
  internal_force_and_stiffness!(domain.assembler, U, domain, Uu, state, props, domain.coords)
  solver.residual .= @views domain.assembler.residuals[domain.dof.unknown_dofs]
  solver.stiffness = SparseArrays.sparse!(domain.assembler)
end

function solve!(solver::LinearSolver)
  solver.unknowns .= zero(eltype(solver.unknowns))
  # cg!(solver.unknowns, solver.stiffness, solver.residual)
  solver.unknowns .= solver.stiffness \ solver.residual
  # ldiv!(solver.unknowns, solver.stiffness, solver.residual)
  # solver.unknowns .= solver.residual
end

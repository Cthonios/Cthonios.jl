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
  K = stiffness(domain, Uu)

  # setup solver
  if :type in keys(input_settings)
    if input_settings[:type] == "default"
      solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
    else
      solver = eval(Meta.parse(input_settings[:type]))()
    end
  else
    solver = LinearSolve.DefaultLinearSolver(LinearSolve.DefaultAlgorithmChoice.UMFPACKFactorization)
  end

  # set up preconditioner
  if :preconditioner in keys(input_settings)
    if input_settings[:preconditioner] == "default"
      Pl = IdentityOperator(length(Uu))
    elseif input_settings[:preconditioner] == "cholesky"
      Pl = cholesky(K)
    else
      Pl = eval(Meta.parse(input_settings[:preconditioner]))(K)
    end
  else
    Pl = IdentityOperator(length(Uu))
  end

  return LinearSolver(Pl, R, K, Uu)
end

function Base.show(io::IO, solver::LinearSolver)
  print(io, "LinearSolver\n", 
            "  Preconditioner = $(solver.preconditioner)")
end

function update_residual!(
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}
  residual!(domain, Uu, domain.coords, U)
  solver.residual .= @views domain.assembler.residuals[domain.dof.unknown_dofs]
end

function update_stiffness!(
  solver::LinearSolver, domain::QuasiStaticDomain, 
  Uu::V2, U::V1
) where {V1 <: NodalField, V2 <: AbstractVector}
  stiffness!(domain, Uu, domain.coords, U)
  solver.stiffness = sparse(domain.assembler)
end

function solve!(solver::LinearSolver)
  solver.unknowns .= zero(eltype(solver.unknowns))
  cg!(solver.unknowns, solver.stiffness, solver.residual)
end
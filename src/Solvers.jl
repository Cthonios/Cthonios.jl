abstract type SolverSettings end
abstract type Solver end

struct NewtonSolverSettings <: SolverSettings
  max_steps::Int
end

struct NewtonSolver{V <: AbstractArray{<:Number, 1}} <: Solver
  settings::NewtonSolverSettings
  ΔUu::V
end

function NewtonSolver(input_settings::D, domain::Domain) where D <: Dict
  ΔUu = zeros(Float64, length(domain.dof.unknown_indices))
  return NewtonSolver(NewtonSolverSettings(20), ΔUu)
end

function solve!(domain::Domain, solver::NewtonSolver)
  # op = QuasiStaticDomainOperator{Float64}(domain)
  dof       = domain.dof
  assembler = domain.assembler

  # extra assembly for now since it check pos def at beginning
  # update_fields!(domain)
  # assemble!(domain)

  for n in 1:solver.settings.max_steps
    update_fields!(domain)
    assemble!(domain)

    K = assembler.K[dof.is_unknown, dof.is_unknown]
    R = assembler.R[dof.is_unknown]
    # UpdatePreconditioner!(P, K, 2)

    # P = CholeskyPreconditioner(K, 2)
    # P = AMGPreconditioner{RugeStuben}(K)
    # P = AMGPreconditioner{SmoothedAggregation}(K)
    # TODO genrealize solver
    # cg!(solver.ΔUu, K, R; Pl=P)
    # gmres!(solver.ΔUu, K, R; Pl=P)
    # update unknowns

    solver.ΔUu .= K \ R

    domain.Uu .= domain.Uu .- solver.ΔUu

    norm_R = norm(domain.assembler.R[domain.dof.is_unknown])
    norm_U = norm(solver.ΔUu)

    # TODO redo
    @info @sprintf "Iteration %5i: ||R|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R norm_U

    # TODO make solver tolerance a setting
    if norm_R < 1e-12 || norm_U < 1e-12
      break
    end
  end

end

include("TrustRegionSolver.jl")
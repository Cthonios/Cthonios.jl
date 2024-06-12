abstract type AbstractLinearSolver end
abstract type AbstractNonlinearSolver end
function setup_solver_cache end

struct DirectSolver <: AbstractLinearSolver
end

function solve(::DirectSolver, obj, Uu)
  R = gradient(obj, Uu)
  K = hessian(obj, Uu)
  ΔUu = -K \ R
  return ΔUu
end

struct SolverCache{Energy, Energies, Asm}
  Π::Energy
  Πs::Energies
  assembler::Asm
end

function Base.similar(cache::SolverCache)
  Π = similar(cache.Π)
  Πs = similar(cache.Πs)
  asm = similar(cache.assembler)
  return SolverCache(Π, Πs, asm)
end

struct NewtonSolver{LinSolver} <: AbstractNonlinearSolver
  linear_solver::LinSolver
end

function logger(::NewtonSolver, n, norm_R, norm_R0, norm_U)
  @info @sprintf "  Iteration %5i: ||R|| = %1.6e    ||R/R0|| = %1.6e    ||ΔUu|| = %1.6e" n norm_R (norm_R / norm_R0) norm_U
end

function solve!(solver::NewtonSolver, obj, Uu)
  tol = 1e-10
  @info "Newton Solver:"
  R0 = norm(residual(obj.domain, Uu))
  for n in 1:10
    ΔUu = solve(solver.linear_solver, obj, Uu)
    Uu .+= ΔUu

    norm_R = norm(obj.domain.cache.solver_cache.assembler.residuals[obj.domain.static.dof.unknown_dofs])
    norm_U = norm(ΔUu)
    logger(solver, n, norm_R, R0, norm_U)
    if norm_R < tol || norm_U < tol
      break
    end
  end
  # return Uu
  return nothing
end

function setup_solver_cache(::NewtonSolver, static)
  assembler = StaticAssembler(static.dof, map(x -> x.fspace, values(static.sections)))
  Πs = Dict{Symbol, Any}()
  for (section_name, section) in pairs(static.sections)
    NQ, NE = FiniteElementContainers.num_q_points(section.fspace), num_elements(section.fspace)
    Πs[section_name] = ElementField{NQ, NE, Vector, Float64}(undef)
    Πs[section_name] .= zero(Float64)
  end
  return SolverCache(
    zeros(Float64, 1), ComponentArray(Πs),
    assembler
  )
end


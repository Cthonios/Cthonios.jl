# hacks for now

struct EigenSolver{O, A, T, F} <: AbstractSolver
  objective::O
  mass_assembler::A
  stiffness_assembler::A
  timer::T
  nev::Int
  tolerance::F
end

function EigenSolver(obj::AbstractObjective, p, timer, nev)
  @timeit timer "EigenSolver - setup" begin
    mass_asm = StaticAssembler(obj.domain)
    stiffness_asm = StaticAssembler(obj.domain)
    update_unknown_dofs!(obj.domain, mass_asm)
    update_unknown_dofs!(obj.domain, stiffness_asm)
    return EigenSolver(obj, mass_asm, stiffness_asm, timer, nev, 1.e-8)
  end
end

function solve!(solver::EigenSolver, obj, Uu, p)
  # TODO the method below isn't defined
  K = hessian!(solver.stiffness_assembler, obj, Uu, p)
  M = mass_matrix!(solver.mass_assembler, obj, Uu, p)
  results = lobpcg(M, K, true, solver.nev; maxiter=1000)
  return results
end

timer(s::EigenSolver) = s.timer

export EigenSolver

module CthoniosAdaptExt

using Adapt
using Cthonios

# boundary conditions
function Adapt.adapt_structure(to, bc::Cthonios.DirichletBCInternal)
  nodes = Adapt.adapt_structure(to, bc.nodes)
  dofs = Adapt.adapt_structure(to, bc.dofs)
  func = Adapt.adapt_structure(to, bc.func)
  return Cthonios.DirichletBCInternal(nodes, dofs, func)
end

# domains
function Adapt.adapt_structure(to, domain::Cthonios.Domain)
  coords = Adapt.adapt_structure(to, domain.coords)
  dof = Adapt.adapt_structure(to, domain.dof)
  sections = Adapt.adapt_structure(to, domain.sections)
  # TODO what do we need to do with the BCs
  dirichlet_bcs = Adapt.adapt_structure(to, domain.dirichlet_bcs)
  dirichlet_dofs = Adapt.adapt_structure(to, domain.dirichlet_dofs)
  return Cthonios.Domain(coords, dof, sections, dirichlet_bcs, dirichlet_dofs)
end

# objectives
function Adapt.adapt_structure(to, obj::Cthonios.Objective)
  domain = Adapt.adapt_structure(to, obj.domain)
  value = Adapt.adapt_structure(to, obj.value)
  gradient = Adapt.adapt_structure(to, obj.gradient)
  hessian = Adapt.adapt_structure(to, obj.hessian)
  timer = Adapt.adapt_structure(to, obj.timer)
  return Cthonios.Objective(domain, value, gradient, hessian, timer)
end

function Adapt.adapt_structure(to, p::Cthonios.ObjectiveParameters)
  X = Adapt.adapt_structure(to, p.X)
  t = Adapt.adapt_structure(to, p.t)
  Ubc = Adapt.adapt_structure(to, p.Ubc)
  props = Adapt.adapt_structure(to, p.props)
  U = Adapt.adapt_structure(to, p.U)
  hvp_scratch = Adapt.adapt_structure(to, p.hvp_scratch)
  q_vals_scratch = Adapt.adapt_structure(to, p.q_vals_scratch)
  return ObjectiveParameters(X, t, Ubc, props, U, hvp_scratch, q_vals_scratch)
end

# problems
function Adapt.adapt_structure(to, prob::Cthonios.QuasiStaticProblem)
  objective = Adapt.adapt_structure(to, prob.objective)
  solver = Adapt.adapt_structure(to, prob.solver)
  pp = prob.post_processor
  timer = Adapt.adapt_structure(to, prob.timer)
  return Cthonios.QuasiStaticProblem(objective, solver, pp, timer)
end

# properties
function Adapt.adapt_structure(to, props::Cthonios.MaterialProperties)
  props = Adapt.adapt_structure(to, props.props)
  return Cthonios.MaterialProperties(props)
end

# sections
function Adapt.adapt_structure(to, sec::Cthonios.SectionInternal)
  # block_name = Adapt.adapt_structure(to, sec.block_name)
  fspace = Adapt.adapt_structure(to, sec.fspace)
  # physics = Adapt.adapt_structure(to, sec.physics)
  physics = sec.physics # nothing here is really on the GPU
  props = Adapt.adapt_structure(to, sec.props)
  return Cthonios.SectionInternal(fspace, physics, props)
end 

# solvers
function Adapt.adapt_structure(to, solver::Cthonios.DirectSolver)
  assembler = Adapt.adapt_structure(to, solver.assembler)
  # TODO there is not adapt ext for LinearSolve.jl
  # might need to specialize based on backend
  linsolve = Adapt.adapt_structure(to, solver.linsolve)
  timer = Adapt.adapt_structure(to, solver.timer)
  return Cthonios.DirectSolver(assembler, linsolve, timer)
end

function Adapt.adapt_structure(to, solver::Cthonios.NewtonSolver)
  linear_solver = Adapt.adapt_structure(to, solver.linear_solver)
  objective = Adapt.adapt_structure(to, solver.objective)
  ΔUu = Adapt.adapt_structure(to, solver.ΔUu)
  timer = Adapt.adapt_structure(to, solver.timer)
  max_iter = Adapt.adapt_structure(to, solver.max_iter)
  abs_tol = Adapt.adapt_structure(to, solver.abs_tol)
  rel_tol = Adapt.adapt_structure(to, solver.rel_tol)
  use_warm_start = Adapt.adapt_structure(to, solver.use_warm_start)
  return NewtonSolver(linear_solver, objective, ΔUu, timer, max_iter, abs_tol, rel_tol, use_warm_start)
end

# time steppers
function Adapt.adapt_structure(to, ts::T) where T <: Cthonios.QuasiStatic
  start_time = Adapt.adapt_structure(to, ts.start_time)
  end_time = Adapt.adapt_structure(to, ts.end_time)
  current_time = Adapt.adapt_structure(to, ts.current_time)
  current_time_step = Adapt.adapt_structure(to, ts.current_time_step)
  Δt = Adapt.adapt_structure(to, ts.Δt)
  return QuasiStatic(start_time, end_time, current_time, current_time_step, Δt)
end

end # module

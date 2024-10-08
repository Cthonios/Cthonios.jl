@testset ExtendedTestSet "Problems - NewtonSolver" begin
  # setup
  func_1(x, t) = -0.125 * t
  func_2(x, t) = 0.0
  bcs = [
    DirichletBC("nset_outer_bottom", [1, 2], func_2)
    DirichletBC("nset_outer_top", [1], func_2)
    DirichletBC("nset_outer_top", [2], func_1)
  ]
  sections = Section[
    Section(
      Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
      "unnamed_block_1", 2
    )
  ]
  domain = Domain("window_pain_tri3.g", sections, bcs, 2)
  asm = Cthonios.StaticAssembler(domain)
  Cthonios.update_unknown_dofs!(domain, asm)
  coords = coordinates(domain.mesh)
  coords = NodalField{size(coords), Vector}(coords)
  time = ConstantTimeStepper(0.0, 1.0, 0.05)

  # constructors
  timer = TimerOutput()
  objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, timer)
  Uu = Cthonios.create_unknowns(objective.domain)
  p = ObjectiveParameters(objective, time)
  solver = NewtonSolver(objective, p, DirectSolver, timer)

  Uu = Cthonios.create_unknowns(solver)

  # pp
  pp = ExodusPostProcessor("window_pain_tri3.g", "output.e", ["displ_x", "displ_y"])

  problem = QuasiStaticProblem(objective, solver, pp, timer)
  Cthonios.solve!(problem, Uu, p)
end

@testset ExtendedTestSet "Problems - TrustRegionSolver" begin
  # setup
  func_1(x, t) = -0.125 * t
  func_2(x, t) = 0.0
  bcs = [
    DirichletBC("nset_outer_bottom", [1, 2], func_2)
    DirichletBC("nset_outer_top", [1], func_2)
    DirichletBC("nset_outer_top", [2], func_1)
  ]
  sections = Section[
    Section(
      Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
      "unnamed_block_1", 2
    )
  ]
  domain = Domain("window_pain_tri3.g", sections, bcs, 2)
  asm = Cthonios.StaticAssembler(domain)
  Cthonios.update_unknown_dofs!(domain, asm)
  coords = coordinates(domain.mesh)
  coords = NodalField{size(coords), Vector}(coords)
  time = ConstantTimeStepper(0.0, 1.0, 0.05)

  # constructors
  timer = TimerOutput()
  objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, timer)
  Uu = Cthonios.create_unknowns(objective.domain)
  p = ObjectiveParameters(objective, time)

  solver = TrustRegionSolver(objective, p, timer; use_warm_start=false)

  Uu = Cthonios.create_unknowns(solver)

  # pp
  pp = ExodusPostProcessor("window_pain_tri3.g", "output.e", ["displ_x", "displ_y"])

  problem = QuasiStaticProblem(objective, solver, pp, timer)
  Cthonios.solve!(problem, Uu, p)
end

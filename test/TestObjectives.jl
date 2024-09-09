@testset ExtendedTestSet "Objectives" begin
  func_1(x, t) = -0.5 * t
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
  times = ConstantTimeStepper(0.0, 1.0, 0.0125)
  objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, TimerOutput())
  p = ObjectiveParameters(objective, times)

end
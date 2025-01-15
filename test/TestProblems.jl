# @testset ExtendedTestSet "Problems - Problem - NewtonSolver" begin
#   # setup
#   func_1(x, t) = -0.125 * t
#   func_2(x, t) = 0.0
#   dbcs = [
#     DirichletBC("nset_outer_bottom", [1, 2], func_2)
#     DirichletBC("nset_outer_top", [1], func_2)
#     DirichletBC("nset_outer_top", [2], func_1)
#   ]
#   nbcs = [

#   ]
#   sections = Section[
#     Section(
#       "unnamed_block_1", 2,
#       Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
#       MaterialProperties(
#         "bulk modulus"  => 0.833,
#         "shear modulus" => 0.3846
#       )
#     )
#   ]
#   domain = Domain("window_pain_tri3.g", sections, dbcs, nbcs)
#   asm = Cthonios.StaticAssembler(domain)
#   Cthonios.update_unknown_dofs!(domain, asm)
#   coords = domain.coords
#   integrator = QuasiStatic(0.0, 1.0, 0.05)

#   # constructors
#   timer = TimerOutput()
#   objective = Objective(domain, Cthonios.energy)
#   pp = ExodusPostProcessor("window_pain_tri3.g", "output.e", ["displ_x", "displ_y"])
#   problem = Problem(objective, integrator, NewtonSolver, pp)
#   Uu, p = Cthonios.create_unknowns_and_parameters(problem)
#   Cthonios.solve!(problem, Uu, p)
# end

@testset ExtendedTestSet "Problems - Problem - TrustRegionSolver" begin
  # setup
  func_1(x, t) = -0.125 * t
  func_2(x, t) = 0.0
  dbcs = [
    DirichletBC("nset_outer_bottom", [1, 2], func_2)
    DirichletBC("nset_outer_top", [1], func_2)
    DirichletBC("nset_outer_top", [2], func_1)
  ]
  nbcs = [

  ]
  sections = Section[
    Section(
      "unnamed_block_1", 2,
      Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
      MaterialProperties(
        "bulk modulus"  => 0.833,
        "shear modulus" => 0.3846
      )
    )
  ]
  domain = Domain("window_pain_tri3.g", sections, dbcs, nbcs)
  asm = Cthonios.StaticAssembler(domain)
  Cthonios.update_unknown_dofs!(domain, asm)
  coords = domain.coords
  integrator = QuasiStatic(0.0, 1.0, 0.05)

  # constructors
  timer = TimerOutput()
  objective = UnconstrainedObjective(domain, Cthonios.energy)
  pp = ExodusPostProcessor("window_pain_tri3.g", "output.e", ["displ_x", "displ_y"])
  problem = Problem(objective, integrator, TrustRegionSolver, pp)
  Uu, p = Cthonios.create_unknowns_and_parameters(problem)
  Cthonios.solve!(problem, Uu, p)
end

# @testset ExtendedTestSet "Problems - EigenProblem - EigenSolver" begin
#   # setup
#   func_1(x, t) = -0.125 * t
#   func_2(x, t) = 0.0
#   dbcs = [
#     DirichletBC("nset_outer_bottom", [1, 2], func_2)
#     DirichletBC("nset_outer_top", [1], func_2)
#     DirichletBC("nset_outer_top", [2], func_1)
#   ]
#   nbcs = [

#   ]
#   sections = Section[
#     Section(
#       "unnamed_block_1", 2,
#       Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
#       MaterialProperties(
#         "bulk modulus"  => 0.833,
#         "shear modulus" => 0.3846
#       )
#     )
#   ]
#   domain = Domain("window_pain_tri3.g", sections, dbcs, nbcs)
#   asm = Cthonios.StaticAssembler(domain)
#   Cthonios.update_unknown_dofs!(domain, asm)
#   coords = domain.coords
#   integrator = QuasiStatic(0.0, 1.0, 0.05)

#   # constructors
#   timer = TimerOutput()
#   objective = Objective(domain, Cthonios.energy)
#   pp = ExodusPostProcessor("window_pain_tri3.g", "output.e", ["displ_x", "displ_y"])
#   problem = Problem(objective, integrator, NewtonSolver, pp)
#   Uu, p = Cthonios.create_unknowns_and_parameters(problem)
#   Cthonios.solve!(problem, Uu, p)

#   # eigen problem
#   sections = Section[
#     Section(
#       "unnamed_block_1", 2,
#       Cthonios.SolidDynamics(NeoHookean(), PlaneStrain()),
#       MaterialProperties(
#         "bulk modulus"  => 0.833,
#         "shear modulus" => 0.3846
#       )
#     )
#   ]
#   domain = Domain("window_pain_tri3.g", sections, dbcs, nbcs)
#   objective = Objective(
#     domain, 
#     Cthonios.energy, Cthonios.gradient, Cthonios.hessian, 
#     Cthonios.neumann_energy, Cthonios.neumann_gradient, Cthonios.neumann_hessian, 
#     timer
#   )  
#   p = ObjectiveParameters(objective, time)
#   solver = EigenSolver(objective, p, timer, 10)
#   pp = ExodusPostProcessor("window_pain_tri3.g", "eigen_output.e", ["displ_x", "displ_y"])
#   problem = EigenProblem(objective, solver, pp, timer)
#   Cthonios.solve!(problem, Uu, p)
# end

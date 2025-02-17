@testset ExtendedTestSet "Objectives" begin
  # setup
  func_1(x, t) = -0.5 * t
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
  time = QuasiStatic(0.0, 1.0, 0.0125)

  # constructors
  objective = UnconstrainedObjective(
    domain, Cthonios.energy
    # Cthonios.energy, Cthonios.gradient, Cthonios.hessian, 
    # Cthonios.neumann_energy, Cthonios.neumann_gradient, Cthonios.neumann_hessian, 
    # TimerOutput()
  )
  Uu = Cthonios.create_unknowns(objective.domain)
  p = ObjectiveParameters(objective, time)

  # some methods
  # Cthonios.step!(p)
  Cthonios.update_dirichlet_vals!(p, objective)

  for bc in domain.dirichlet_bcs
    vals = map(x -> bc.func(x, time.current_time[1]), coords[:, bc.nodes])
    for val in vals
      @test val in p.Ubc
    end
  end

  # just testing to see if they run for now
  val = zeros(1)
  g = Cthonios.create_fields(domain)
  Hv = Cthonios.create_fields(domain)
  Vv = Cthonios.create_unknowns(domain)
  Cthonios.objective!(val, objective, Uu, p)
  Cthonios.gradient!(g, objective, Uu, p)
  Cthonios.hvp!(Hv, objective, Uu, p, Vv)
  Cthonios.hessian!(asm, objective, Uu, p)

end
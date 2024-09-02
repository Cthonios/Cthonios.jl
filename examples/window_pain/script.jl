using ConstitutiveModels
using Cthonios
using Exodus
using FiniteElementContainers
using LinearAlgebra

# file management
mesh_file = Base.source_dir() * "/window_pain_tri3.g"

# function main()
# global setup
times = ConstantTimeStepper(0.0, 1.0, 0.0125)
n_dofs = 2

# functions
func_1(x, t) = -0.5 * t
func_2(x, t) = 0.0

# bcs
disp_bcs = [
  DirichletBC("nset_outer_bottom", [1, 2], func_2)
  DirichletBC("nset_outer_top", [1], func_2)
  DirichletBC("nset_outer_top", [2], func_1)
]
traction_bcs = [
]

# sections
sections = TotalLagrangeSection[
  TotalLagrangeSection(
    Cthonios.SolidMechanics(NeoHookean(), PlaneStrain()),
    "unnamed_block_1", 2
  )
]
domain = Domain(mesh_file, times, sections, disp_bcs, 2)
objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian)
solver = NewtonSolver(objective, DirectSolver)

Uu = Cthonios.create_unknowns(solver)
p = nothing

# pp
U = Cthonios.create_fields(domain)
pp = Cthonios.ExodusPostProcessor(mesh_file, "output.e", ["displ_x", "displ_y"])

try 
  # pp
  Cthonios.update_bcs!(U, solver.objective.domain)
  Cthonios.update_fields!(U, domain, Uu)
  Cthonios.write_time(pp, 1, 0.0)
  Cthonios.write_fields(pp, U, 1)

  # loop over steps
  for n in 1:80
    # load step
    Cthonios.step!(solver.objective.domain)
    Cthonios.update_bcs!(solver.linear_solver.U, solver.objective.domain)
    Cthonios.solve!(solver, Uu, p)

    # pp
    Cthonios.update_bcs!(U, solver.objective.domain)
    Cthonios.update_fields!(U, domain, Uu)
    Cthonios.write_time(pp, n + 1, solver.objective.domain.t.current_time)
    Cthonios.write_fields(pp, U, n + 1)
  end
catch e
  Cthonios.close(pp)
  throw(e)
end

Cthonios.close(pp)

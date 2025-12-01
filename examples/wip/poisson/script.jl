using Cthonios
using Exodus
using IterativeSolvers
using LinearAlgebra
using Parameters
using SparseArrays
using TimerOutputs

f(X, _) = 2. * π^2 * sin(2 * π * X[1]) * sin(2 * π * X[2])
zero_func(_, _) = 0.0

# set up initial containers
mesh_file = Base.source_dir() * "/poisson.g"
output_file = Base.source_dir() * "/output.e"
timer = TimerOutput()

times = ConstantTimeStepper(0.0, 0.0, 0.0 / 1.)
dbcs = DirichletBC[
  DirichletBC("nset_1", [1], zero_func)
  DirichletBC("nset_2", [1], zero_func)
  DirichletBC("nset_3", [1], zero_func)
  DirichletBC("nset_4", [1], zero_func)
]
nbcs = [

]
sections = [
  Section("block_1", 2, Poisson(f), MaterialProperties())
]
domain = Domain(mesh_file, sections, dbcs, nbcs, 1)
objective = Objective(domain, Cthonios.energy, Cthonios.gradient, Cthonios.hessian, timer)

p = ObjectiveParameters(objective, times)

solver = NewtonSolver(objective, p, DirectSolver, timer)
Uu = Cthonios.create_unknowns(solver)
@time Cthonios.solve!(solver, Uu, p)

# # postprocess
# U = Cthonios.create_fields(domain)
# Cthonios.update_fields!(U, domain, Uu)
# prob = Problem(solver) # add post-processor and other top level stuff here

copy_mesh(mesh_file, output_file)
exo = ExodusDatabase(output_file, "rw")
write_names(exo, NodalVariable, ["u"])
write_time(exo, 1, 0.0)
write_values(exo, NodalVariable, 1, "u", p.U[1, :])
close(exo)

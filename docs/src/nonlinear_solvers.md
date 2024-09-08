```@meta
CurrentModule = Cthonios
```
# Nonlinear Solvers
One of our goals is to make adding new and hacking existing nonlinear solvers easy!

# Input File Syntax Example
Below is an example input file block for a trust region solver that leverages a direct solver with an LDLT factorization. The nonlinear solver is also utilizing warm start to improve the first few iterations of each load step.
```
linear solvers:
  direct:
    type: DirectLinearSolver
    factorization method: ldl

nonlinear solvers:
  trs:
    type: TrustRegionSolver
    linear solver: direct
    warm start: on
```

# General methods and abstract types
```@autodocs
Modules = [Cthonios]
Pages = ["solvers/nonlinear_solvers/NonlinearSolvers.jl"]
Order = [:type, :function]
```

# Newton Solver
```@autodocs
Modules = [Cthonios]
Pages = ["solvers/nonlinear_solvers/NewtonSolver.jl"]
Order = [:type, :function]
```

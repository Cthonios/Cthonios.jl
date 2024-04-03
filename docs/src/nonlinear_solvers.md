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

# Newton Solver
```@autodocs
Modules = [Cthonios]
Pages = ["solvers/NewtonSolver.jl"]
Order = [:type, :function]
```

# Trust Region Solver
```@autodocs
Modules = [Cthonios]
Pages = ["solvers/TrustRegionSolver.jl"]
Order = [:type, :function]
```

# General
```@autodocs
Modules = [Cthonios]
Pages = ["solvers/NonlinearSolvers.jl"]
Order = [:type, :function]
```
```@meta
CurrentModule = Cthonios
```
# Linear Solvers
This section describes the interface for linear solvers.

# Input File Syntax Example
Below is an example input file block for a direct solver with an LDLT factorization.
```
linear solvers:
  direct:
    type: DirectLinearSolver
    factorization method: ldl
```

# Direct Solver
```@autodocs
Modules = [Cthonios]
Pages = ["solvers/linear_solvers/DirectLinearSolver.jl"]
Order = [:type, :function]
```

# General
```@autodocs
Modules = [Cthonios]
Pages = ["solvers/linear_solvers/LinearSolvers.jl"]
Order = [:type, :function]
```

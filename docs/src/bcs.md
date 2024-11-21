# Boundary Conditions

## Abstract Types
There is currently only a loose set of abstract types to define BCs.
TODO
```@autodocs
Modules = [Cthonios]
Pages = ["bcs/BoundaryConditions.jl"]
Order = [:type, :function]
```

## Dirichlet BCs
Dirichlet BCs are defined in a general way where each ```DirichletBC``` acts on a collection of nodes and a collection of dof indices.

We can set up a Dirichlet BC in a script with syntax similar to the following
```jldoctest
using Cthonios
bc = DirichletBC("nset", [1], (x, t) -> 0.0)

# output
DirichletBC:
  Node set = nset
  Dofs     = [1]
  Function = #1
```
In general, problems will have a collection of Dirichlet BCs which can be set up as follows
```jldoctest
using Cthonios
dbcs = [
  DirichletBC("nset_1", [1, 2], (x, t) -> 0.0)
  DirichletBC("nset_2", [1], (x, t) -> 0.0)
  DirichletBC("nset_2", [2], (x, t) -> 1.0)
]

# output
3-element Vector{DirichletBC{String, Vector{Int64}}}:
 DirichletBC:
  Node set = nset_1
  Dofs     = [1, 2]
  Function = #1

 DirichletBC:
  Node set = nset_2
  Dofs     = [1]
  Function = #2

 DirichletBC:
  Node set = nset_2
  Dofs     = [2]
  Function = #3

```

Internally these types defined in a script will be converted to a new type ```DirichletBCInternal``` which will read quantities off of the input mesh. This allows for BC definition prior to mesh reading so error checking can be done more efficiently.

```@autodocs
Modules = [Cthonios]
Pages = ["bcs/DirichletBC.jl"]
Order = [:type, :function]
```

## Neumann BCs
Neumann BCs are nearly identical to Dirichlet BCs in terms of script inputs except these BCs do not require a set of dofs. One can define a NeumannBC as follows

```jldoctest
using Cthonios
bc = NeumannBC("sset", (x, t) -> @SVector [0., 1.])

# output
NeumannBC:
  Side set = sset
  Function = #1
```

Collections of ```NeumannBC```s follows identically to collections of ```DirichletBC```s.

```@autodocs
Modules = [Cthonios]
Pages = ["bcs/NeumannBC.jl"]
Order = [:type, :function]
```

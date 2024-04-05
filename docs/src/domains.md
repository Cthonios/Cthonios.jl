```@meta
CurrentModule = Cthonios
```
# Domains
A domain in Cthonios is a group of element sections and all the associated datastructures necessary to run a finite element calculation. Different types of domains are used for different problems, e.g. quasi-static vs. dynamics problems have different requirements and solution strategies.

# Input file syntax example
```
domains:
  domain_1:
    mesh:
      ...
    boundary conditions:
      ...
    sections:
      ...
    time stepper:
      ...
```

# Quasi-static Domain
Quasi-static domains are appropriate for problems where the following equation holds

``
\nabla\cdot \mathbf{P} + \mathbf{b} = \mathbf{0},
``
where ``\mathbf{P}`` is the first Piola-Kirchoff stress, ``\nabla`` is the gradient operator with respsect to the reference corrdinates, and ``\mathbf{b}`` is a body force.

# Index
```@autodocs
Modules = [Cthonios]
Pages = ["domains/QuasiStaticDomain.jl"]
Order = [:type, :function]
```
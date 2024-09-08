```@meta
CurrentModule = Cthonios
```
# Sections
Within Cthonios, sections refer to groups or "blocks" of elements that share a common element topology and material model. We allow by default for each element to possess different fixed material properties and an arbitrary number of state variables. Since the material model is the same for all elements in a section, the number of state variables will be the same at each element/quadrature point.

# Input file syntax
TODO

```@autodocs
Modules = [Cthonios]
Pages = ["Sections.jl"]
Order = [:type, :function]
```

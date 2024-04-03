```@meta
CurrentModule = Cthonios
```
# Materials
Cthonios utilizes [```ConstitutiveModels.jl```](https://github.com/Cthonios/ConstitutiveModels.jl) as a material model backend.

Here is an example for how to define multiple materials within an input file
```
materials:
  metal_linear:
    model: LinearElastic
    properties:
      youngs modulus: 1.0
      poissons ratio: 0.3
  soft rubber:
    model: NeoHookean
    properties:
      bulk modulus: 50.0
      shear modulus: 1.0
```

# Useful methods
```@docs
setup_material
```

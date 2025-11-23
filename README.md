# Cthonios

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cthonios.github.io/Cthonios.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cthonios.github.io/Cthonios.jl/dev/)
[![Build Status](https://github.com/Cthonios/Cthonios.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Cthonios/Cthonios.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Cthonios/Cthonios.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Cthonios/Cthonios.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/C/Cthonios.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/C/Cthonios.html)

Cthonios - Continuum THermodynamic Optimization-based Numerical Implementations Of Solid mechanics

This is a collection of tools for solving solid mechanics problems via optimization methods. This is larrgely inspired by the python package [OptimiSM](https://github.com/sandialabs/optimism).

The tools build on several other packages such as [Exodus.jl](https://github.com/cthonios/Exodus.jl), [ReferenceFiniteElements.jl](https://github.com/Cthonios/ReferenceFiniteElements.jl), [FiniteElementContainers.jl](https://github.com/Cthonios/FiniteElementContainers.jl), and [ConstitutiveModels.jl](https://github.com/Cthonios/ConstitutiveModels.jl) to name a few of the direct dependencies developed by the [Cthonios organization](https://github.com/Cthonios).

# Table of Contents
1. [Installation](#installation)
2. [Documentation](#documentation)
3. [Usage as an application](#usage-as-an-application)
4. [Example Script Solving a Forward Problem](#example-script-solving-a-forward-problem)
5. [Example Input File Solving a Forward Problem](#example-input-file-solving-a-forward-problem)
6. [Contributing](#contributing)

## Installation
From the julia package manager (which can be entered by typing the "]" key after starting up julia) one can do the following
```
pkg> add Cthonios
```

Or one can use the ```Pkg``` module inside julia as follows
```julia
using Pkg
Pkg.add("Cthonios")
```

## Documentation
Documentation deployed from main branch can be found [here](https://cthonios.github.io/Cthonios.jl/dev/)

## Usage as an application
Using ```cthonios``` as an application is still very much in development. The easiest route is to do the following

```julia
using Cthonios
push!(ARGS, "-i")
push!(ARGS, "my_input_file.yaml")
Cthonios.cthonios_main(ARGS)
```

To avoid the above, an executable can be built with either [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) or [JuliaC](https://github.com/JuliaLang/JuliaC.jl)(requires julia 1.12)

## Example Script Solving a Forward Problem
Below is an example julia script for solving a column buckling problem in 3D with a trust region solver

```julia
# Load up necessary packages
using ConstitutiveModels
using Cthonios

# File management
mesh_file = Base.source_dir() * "/mesh.g"
output_file = splitext(mesh_file)[1] * "-output.exo"

# Times
times = TimeStepper(0., 1., 40)

# Physics and properties
physics = (;
    block_1 = SolidMechanics(
        ThreeDimensional(), NeoHookean()
    )
)
props = (;
    block_1 = Dict{String, Any}(
      "bulk modulus"  => 10.0,
      "shear modulus" => 1.0
    )
)

# Functions for BCs
func_1(x, t) = 0.0
func_2(x, t) = -0.5 * t
func_3(x, t) = @SVector [0., -0.025 * t, 0.]

# Boundary conditions
dirichlet_bcs = [
  DirichletBC("displ_x", "sset_y_negative", func_1)
  DirichletBC("displ_y", "sset_y_negative", func_1)
  DirichletBC("displ_z", "sset_y_negative", func_1)
  DirichletBC("displ_x", "sset_y_positive", func_1)
  DirichletBC("displ_z", "sset_y_positive", func_1)
  DirichletBC("displ_y", "sset_y_positive", func_2)
]

# Simulation setup
sim = SingleDomainSimulation(
    mesh_file, output_file, 
    times, physics, props;
    dirichlet_bcs=dirichlet_bcs
)

# objective setup
objective = Cthonios.QuasiStaticObjective()
objective_cache = Cthonios.setup_cache(objective, sim)

# solver setup
solver = Cthonios.TrustRegionSolverGPU(objective_cache; use_warm_start=true)

Cthonios.run!(objective_cache, solver, sim)
```

![](./assets/column-buckle-3d.gif)

# Example Input File Solving a Forward Problem
Below is an example input file (currently only YAML is supported) for solving a forward problem of metamaterial structure in 2D.

```yaml
simulation:
  type: SingleDomainSimulation
  forward problem objective: QuasiStaticObjective

  mesh:
    type: UnstructuredMesh
    file name: examples/hole_array/mesh/hole_array_tri6.exo

  time:
    start time: 0.0
    end time: 1.0
    number of steps: 80

  materials:
    mat_1:
      NeoHookean:
        density: 1.0
        Young's modulus: 1.
        Poisson's ratio: 0.495
      Gent:
        density: 1.0
        Young's modulus: 1.
        Poisson's ratio: 0.495
        Jm: 3.

  physics:
    type: SolidMechanics
    kinematics formulation: PlaneStrain
    material assignment:
      mat_1:
        blocks: [Block1]
        model: NeoHookean

  nonlinear solvers:
    trs:
      type: TrustRegionSolverGPU
      verbose: true
      use_warm_start: true

  initial conditions: []

  dirichlet boundary conditions:
    - fields: ["displ_x", "displ_y"]
      function: (x, t) -> 0.0
      sidesets: [yminus_sideset]
    - fields: ["displ_x"]
      function: (x, t) -> 0.0
      sidesets: [yplus_sideset]
    - fields: ["displ_y"]
      function: (x, t) -> -8. * t
      sidesets: [yplus_sideset]

  neumann boundary conditions: []

  contact pairs: []

  nonlinear solver: trs

```
![](./assets/hole-array-quasi-static.gif)

# Contributing
Please file issues, open PRs, open discussions, etc.
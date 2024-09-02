```@meta
CurrentModule = Cthonios
```
# Running from the REPL
To run ```Cthonios``` from the REPL one can utilize the following method

TODO

# Running Cthonios as an executable
To run ```Cthonios``` as an executable, one must first run the ```CthoniosBuild.jl``` script (assuming one has ```PackageCompiler``` installed already). This will build an executable called ```cthonios``` in a ```build``` folder. To run ```Cthonios``` you can use the following command

```
/path/to/cthonios -i <input-file.yaml>
```

# Complete example input file
```
# Global scope
functions:
  displacement_ramp: 
    type: ScalarFunction{2, Float64, Float64, Float64}
    expression: (x, t) -> -0.5 * t
  zero_func:
    type: ScalarFunction{2, Float64, Float64, Float64}
    expression: (x, t) -> 0.0

materials:
  soft rubber:
    model: NeoHookean
    properties:
      bulk modulus: 50.0
      shear modulus: 1.0

linear solvers:
  direct:
    type: DirectLinearSolver
    factorization method: ldl

nonlinear solvers:
  trs:
    type: TrustRegionSolver
    linear solver: direct
    warm start: on

# Domains scope
domains:
  domain_1:
    mesh:
      type: ExodusDatabase{Int32, Int32, Int32, Float64}
      file name: window_pain_tri3.g

    boundary conditions:
      displacement: 
      - nodeset ids: [3]
        dofs: [1, 2]
        function: zero_func
      - nodeset ids: [1]
        dofs: [1]
        function: zero_func
      - nodeset ids: [1]
        dofs: [2]
        function: displacement_ramp

    sections:
    - type: TotalLagrangeSection
      block id: 1
      formulation: plane strain
      material: soft rubber

    time stepper:
      type: ConstantTimeStepper
      start time: 0.0
      end time: 1.0
      time step: 0.025

problems:
  - type: ForwardProblem
    domain: domain_1
    solver: trs
    results:
      output file name: output.e
      nodal fields:
      - displacement
      - internal force
      element fields:
      - properties
      quadrature fields:
      - state variables

  - type: EnergySensitivityProblem
    domain: domain_1
    solver: trs
    results:
      output file name: gradients.e
      nodal fields:
      - displacement
      - dcoordinates
      element fields:
      - dproperties
```
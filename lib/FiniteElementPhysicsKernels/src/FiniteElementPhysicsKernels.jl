module FiniteElementPhysicsKernels

export setup_function
export stable_time_step
export HeatConduction
export SolidMechanics

import FiniteElementContainers
import FiniteElementContainers: AbstractPhysics,
    FunctionSpace,
    GeneralFormulation,
    discrete_gradient,
    extract_stiffness,
    extract_stress,
    interpolate_field_gradients,
    interpolate_field_values_and_gradients,
    interpolate_field_values,
    map_interpolants,
    modify_field_gradients,
    num_dimensions,
    num_fields,
    scatter_with_values!,
    scatter_with_gradients!,
    scatter_with_gradients_and_gradients!,
    scatter_with_values_and_values!
using ConstitutiveModels
using StaticArrays
using Tensors

include("Common.jl")

# physics implementations
include("HeatConduction.jl")
include("LinearSolidMechanics.jl")
include("SolidMechanics.jl")

end # module FiniteElementPhysicsKernels

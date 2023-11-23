module Cthonios

# edpxort
export CthoniosCommon
export DofManager
export Domain
export EssentialBC


# dev stuff
export Connectivity
export connectivity
export dof_connectivity

export NonAllocatedFunctionSpace
export PreAllocatedFunctionSpace
export quadrature_point
export shape_function_gradients
export shape_function_values

export TotalLagrangeSection

# dependencies
using ConstitutiveModels
using DocStringExtensions
using Exodus
using LinearAlgebra
using Logging
using LoggingExtras
using PaddedViews
using ReferenceFiniteElements
using StaticArrays
using StructArrays
using Tensors
using TimerOutputs
using YAML

# for docs
@template (FUNCTIONS, METHODS, MACROS) = 
"""
$(TYPEDSIGNATURES)
$(DOCSTRING)
$(METHODLIST)
"""

@template (TYPES) = 
"""
$(TYPEDFIELDS)
$(DOCSTRING)
"""

# high level stuff
include("Common.jl")

# low level containers
include("BoundaryConditions.jl")
include("DofManagers.jl")
include("FunctionSpaces.jl")
include("Sections.jl")

# high level containers
include("Domains.jl")
include("Mechanics.jl")
# include("Connectivities.jl")

# main function, eventually add a CLI wrapper
function cthonios_main(input_file::String)
  log_file_name = splitext(input_file)[1] * ".log"

  common = CthoniosCommon(log_file_name)

  cthonios_header(common)
  dump_input_file(common, input_file)

  domain = Domain(common, input_file)

  # new_section(common, "Timings")
  with_logger(common) do
    new_section("Timings")
    @info timer(common)
  end
end

end # module
module Cthonios

using Aprepro_jll
using ArgParse
using Exodus
using FiniteElementContainers
using LinearAlgebra
using Logging
using LoggingExtras
using ReferenceFiniteElements
using SparseArrays
using Suppressor
using YAML

# upper level abstract types
abstract type CthoniosException <: Exception end

# lego pieces to build together
include("Aprepro.jl")
include("Headers.jl")
include("Parsers.jl")
include("Systems.jl")

# where the main function is
include("CLI.jl")

end # module
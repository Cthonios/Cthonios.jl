"""
$(TYPEDEF)
"""
abstract type AbstractBCInput end
"""
$(TYPEDEF)
"""
abstract type AbstractBCInternal end

include("DirichletBC.jl")
include("NeumannBC.jl")
include("SchwarzBC.jl")

# exports
export DirichletBC
export NeumannBC

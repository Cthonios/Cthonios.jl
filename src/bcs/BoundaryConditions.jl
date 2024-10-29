abstract type AbstractBCInput end
abstract type AbstractBCInternal end

include("DirichletBC.jl")
include("NeumannBC.jl")

# exports
export DirichletBC
export NeumannBC

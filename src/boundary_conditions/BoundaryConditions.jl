abstract type AbstractBC end
abstract type AbstractDirichletBC{Nodes, Dofs} <: AbstractBC end

include("DisplacementBC.jl")

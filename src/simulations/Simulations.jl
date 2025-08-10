abstract type AbstractSimulation end

abstract type AbstractSimulationCache{
    A  <: FiniteElementContainers.AbstractAssembler, 
    P1 <: Parameters,
    P2 <: FiniteElementContainers.PostProcessor,
    T  <: TimerOutput
} end

function FiniteElementContainers.create_unknowns(sim::AbstractSimulationCache)
    return create_unknowns(sim.assembler, H1Field)
end

function parameters(sim::AbstractSimulationCache)
    return sim.parameters
end

include("SingleDomainSimulation.jl")

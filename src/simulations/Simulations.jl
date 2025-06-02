abstract type AbstractSimulation{
    A <: FiniteElementContainers.AbstractAssembler, 
    P <: Parameters,
    T <: TimerOutput
} end

function FiniteElementContainers.create_unknowns(sim::AbstractSimulation)
    return create_unknowns(sim.assembler, H1Field)
end

function parameters(sim::AbstractSimulation)
    return sim.parameters
end

include("SingleDomainSimulation.jl")

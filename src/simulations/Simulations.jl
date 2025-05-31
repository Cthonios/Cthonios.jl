abstract type AbstractSimulation{
    I <: FiniteElementContainers.AbstractIntegrator, 
    P <: PostProcessor, 
    S <: FiniteElementContainers.AbstractNonLinearSolver, 
    T <: TimerOutput
} end

include("SingleDomainSimulation.jl")

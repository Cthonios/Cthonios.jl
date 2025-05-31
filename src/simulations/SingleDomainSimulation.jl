struct SingleDomainSimulation{I, P, S, T} <: AbstractSimulation{I, P, S, T}
    integrator::I
    post_processor::P
    solver::S
    timer::T
end

function SingleDomainSimulation(
    physics::AbstractPhysics, 
    integrator::AbstractIntegrator, 
    solver::FiniteElementContainers.AbstractNonLinearSolver
)

end

function SingleDomainSimulation(input_settings::Dict{Symbol, Any})

end
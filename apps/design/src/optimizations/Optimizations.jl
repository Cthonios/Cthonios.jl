# TODO more to add....
# keeping lightweight for now
struct SimulationTemplate{
    O,
    T,
    Ph,
    Pr,
    D
}
    objective_type::O
    mesh_file::String
    output_file::String
    times::T
    physics::Ph
    properties::Pr
    dirichlet_bcs::D
end

function run(sim::SimulationTemplate, solver)
    sim = SingleDomainSimulation(
        sim.objective_type,
        sim.mesh_file, sim.output_file,
        sim.times, sim.physics, sim.properties;
        dirichlet_bcs = sim.dirichlet_bcs
    )
    objective = Cthonios.get_objective(sim)
    u, p = Cthonios.get_unknowns_and_parameters(sim)
    solver = solver(objective, u, p)
    Cthonios.run!(solver, objective, u, p, sim.mesh, sim.output_file)
    return sim
end

abstract type AbstractDesignOptimization end
include("StructuralOptimizations.jl")

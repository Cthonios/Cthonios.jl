struct UnconstrainedObjective{A, F1, F2, F3, T} <: AbstractObjectiveCache{A, F1, T}
    assembler::A
    value::F1
    gradient_u::F2
    hessian_u::F3
    timer::T
end

function UnconstrainedObjective(
    sim, value, gradient, hessian
)
    return UnconstrainedObjective(
        sim.assembler, 
        value, gradient, hessian,
        TimerOutput()
    )
end

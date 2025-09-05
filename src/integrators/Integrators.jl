abstract type AbstractTimeIntegrator{
    O, # objective cache
    S, # solver
    T <: TimeStepper
} end

function step! end

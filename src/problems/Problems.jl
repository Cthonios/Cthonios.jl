"""
$(TYPEDEF)
Abstract base problem type.
```O``` - ```Objective````
```S``` - ```Solver```
```P``` - ```PostProcessor```
```T``` - ```TimerOutput```
"""
abstract type AbstractProblem{O, S, P, T} end

"""
$(TYPEDSIGNATURES)
"""
timer(prob::AbstractProblem) = prob.timer

# implementations
include("EigenProblem.jl")
include("QuasiStaticProblem.jl")

# exports
export EigenProblem
export QuasiStaticProblem

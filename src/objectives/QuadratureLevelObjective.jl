struct QuadratureLevelObjective{F1, F2, F3} <: AbstractObjective{F1}
    value::F1
    gradient_u::F2
    hessian_u::F3
end

# struct QuadratureLevelObjectiveCache{O, S, T} <: AbstractObjectiveCache{O, S, T}
#     objective::O
#     sim_cache::S
#     timer::T
# end

# function QuadratureLevelObjectiveCache(
#     objective::QuadratureLevelObjective,
#     sim
# )
#     sim_cache = simulation_cache(sim)
#     return QuadratureLevelObjectiveCache(objective, sim_cache, TimerOutput())
# end

# function gradient(o::QuadratureLevelObjectiveCache, Uu, p)
#     @timeit o.timer "Objective - gradient!" begin
#         # assemble_vector!(assembler(o), Uu, p, H1Field, o.objective.gradient_u)
#         assemble_vector!(assembler(o), o.objective.gradient_u, Uu, p)
#     end
#     return residual(assembler(o))
# end

# function hessian(o::QuadratureLevelObjectiveCache, Uu, p)
#     @timeit o.timer "Objective - hessian!" begin
#         # assemble_matrix!(assembler(o), Uu, p, H1Field, o.objective.hessian_u)
#         assemble_stiffness!(assembler(o), o.objective.hessian_u, Uu, p)
#     end
#     return stiffness(assembler(o)) |> Symmetric # needed for cholesky
# end

# function hvp(o::QuadratureLevelObjectiveCache, Uu, p, Vu)
#     @timeit o.timer "Objective - gradient!" begin
#         # this first one allocates for some reason
#         # assemble_matrix_action!(o.assembler, Uu, p, Vu, H1Field, o.hessian_u)
#         # this one does not
#         # assemble!(assembler(o), Uu, p, Vu, Val{:stiffness_action}(), H1Field)
#         assemble_matrix_action!(assembler(o), o.objective.hessian_u, Uu, Vu, p)
#     end
#     return FiniteElementContainers.hvp(assembler(o))
# end

# function value(o::QuadratureLevelObjectiveCache, Uu, p)
#     @timeit o.timer "Objective - value!" begin
#         assemble_scalar!(assembler(o), o.objective.value, Uu, p)
#     end
#     val = mapreduce(x -> sum(x), sum, assembler(o).scalar_quadarature_storage)
#     return val
# end

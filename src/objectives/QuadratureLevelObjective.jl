struct QuadratureLevelObjective{F1, F2, F3} <: AbstractObjective{F1}
    value::F1
    gradient_u::F2
    hessian_u::F3
end

struct QuadratureLevelObjectiveCache{O, S, T} <: AbstractObjectiveCache{O, S, T}
    # assembler::A
    # cache::C
    objective::O
    sim_cache::S
    timer::T
end

function QuadratureLevelObjectiveCache(
    objective::QuadratureLevelObjective,
    sim
)
    sim_cache = simulation_cache(sim)
    return QuadratureLevelObjectiveCache(objective, sim_cache, TimerOutput())
end

function gradient(o::QuadratureLevelObjectiveCache, Uu, p)
    @timeit o.timer "Objective - gradient!" begin
        assemble_vector!(assembler(o), Uu, p, H1Field, o.objective.gradient_u)
    end
    return residual(assembler(o))
end

function hessian(o::QuadratureLevelObjectiveCache, Uu, p)
    @timeit o.timer "Objective - gradient!" begin
        assemble_matrix!(assembler(o), Uu, p, H1Field, o.objective.hessian_u)
    end
    return stiffness(assembler(o)) |> Symmetric # needed for cholesky
end

function hvp(o::QuadratureLevelObjectiveCache, Uu, p, Vu)
    @timeit o.timer "Objective - gradient!" begin
        # this first one allocates for some reason
        # assemble_matrix_action!(o.assembler, Uu, p, Vu, H1Field, o.hessian_u)
        # this one does not
        assemble!(assembler(o), Uu, p, Vu, Val{:stiffness_action}(), H1Field)
    end
    return FiniteElementContainers.hvp(assembler(o))
end

function value(o::QuadratureLevelObjectiveCache, Uu, p)
    @timeit o.timer "Objective - value!" begin
        assemble_scalar!(assembler(o), Uu, p, H1Field, o.objective.value)
    end
    val = mapreduce(x -> sum(x), sum, assembler(o).scalar_quadarature_storage)
    return val
end

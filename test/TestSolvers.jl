struct MyDummyObjectiveCache1
end

function Cthonios.gradient(::MyDummyObjectiveCache1, x, p)
    return [x[1]^3 - 2x[1] - 5]
end

function Cthonios.hessian(::MyDummyObjectiveCache1, x, p)
    return [3x[1]^2 - x[1]]
end

function test_newton_solver_func_1(verbose)
    objective_cache = MyDummyObjectiveCache1()
    solver = Cthonios.NewtonSolver(objective_cache; verbose=verbose)
    x = ones(1)
    p = nothing
    Cthonios.solve!(solver, x, p)
    @test x[1] ≈ 2.094551481
end

function test_newton_solver_func_1_bad_guess()
    objective_cache = MyDummyObjectiveCache1()
    solver = Cthonios.NewtonSolver(objective_cache; verbose=false)
    x = zeros(1)
    p = nothing
    @test_throws ErrorException Cthonios.solve!(solver, x, p)
    # @test x[1] ≈ 2.094551481
end

struct MyDummyObjectiveCache2
end

function Cthonios.gradient(::MyDummyObjectiveCache2, x, p)
    return [x[1]^2 + x[2]^2 - 4, x[1] - x[2]]
end

function Cthonios.hessian(::MyDummyObjectiveCache2, x, p)
    return [
        2x[1] 2x[2];
        1. -1.
    ]
end

function test_newton_solver_func_2()
    objective_cache = MyDummyObjectiveCache2()
    solver = Cthonios.NewtonSolver(objective_cache)
    x = ones(2)
    p = nothing
    Cthonios.solve!(solver, x, p)
    @test x[1] ≈ sqrt(2)
    @test x[2] ≈ sqrt(2)
end

function test_newton_solver()
    test_newton_solver_func_1(false)
    test_newton_solver_func_1(true)
    test_newton_solver_func_1_bad_guess()
    test_newton_solver_func_2()
end

struct MyDummyObjectiveCache3
end

function Cthonios.value(::MyDummyObjectiveCache3, x, p)
    p = 1.
    return x[1] * (x[1] + 1) + 
           0.3 * x[2] * (x[2] - 0.2) + 
           0.2 * x[3] * (x[3] - 0.5) + 
           x[1] * x[1] * x[2] * x[2] + 
           p * x[1] * x[2] + 
           sin(x[1])

end

function Cthonios.gradient(o::MyDummyObjectiveCache3, x, p)
    return ForwardDiff.gradient(z -> Cthonios.value(o, z, p), x)
end

function Cthonios.hessian(o::MyDummyObjectiveCache3, x, p)
    return ForwardDiff.hessian(z -> Cthonios.value(o, z, p), x)
end

function Cthonios.hvp(o::MyDummyObjectiveCache3, x, v, p)
    return Cthonios.hessian(o, x, p) * v
end

function FiniteElementContainers.create_unknowns(::MyDummyObjectiveCache3)
    return zeros(3)
end

function test_trust_region_solver(verbose)
    objective_cache = MyDummyObjectiveCache3()
    solver = Cthonios.TrustRegionSolver(
        objective_cache;
        preconditioner=Cthonios.NoPreconditioner,
        verbose=verbose
    )
    x = [2., 7., -1.]
    p = nothing
    Cthonios.solve!(solver, x, p)
    g = Cthonios.gradient(objective_cache, x, p)
    @test isapprox(norm(g), 0.; rtol=1.e7)
end

@testset "Newton solver" begin
    test_newton_solver()
end

@testset "Trust region solver" begin
    test_trust_region_solver(false)
    test_trust_region_solver(true)
end

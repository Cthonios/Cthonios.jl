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

@testset "Newton solver" begin
    test_newton_solver()
end

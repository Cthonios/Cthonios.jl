abstract type TestObjective <: Cthonios.AbstractObjective{Float64, Vector{Float64}, Nothing} end

# struct MyDummyObjective1 <: TestObjective
# end

# function Cthonios.gradient(::MyDummyObjective1, x, p)
#     return [x[1]^3 - 2x[1] - 5]
# end

# function Cthonios.hessian(::MyDummyObjective1, x, p)
#     return [3x[1]^2 - x[1]]
# end

# function test_newton_solver_func_1(verbose)
#     objective = MyDummyObjective1()
#     x = ones(1)
#     p = nothing
#     P = NoPreconditioner(objective, x, p)
#     solver = Cthonios.NewtonSolver(objective, x, p; verbose = verbose)
#     Cthonios.solve!(solver, objective, x, p)
#     @test x[1] ≈ 2.094551481
# end

# function test_newton_solver_func_1_bad_guess()
#     objective_cache = MyDummyObjective1()
#     x = zeros(1)
#     p = nothing
#     P = NoPreconditioner(objective, x, p)
#     solver = Cthonios.NewtonSolver(objective, x, p; verbose = false)
#     @test_throws ErrorException Cthonios.solve!(solver, x, p, P)
#     # @test x[1] ≈ 2.094551481
# end

struct MyDummyObjective2 <: TestObjective
end

function Cthonios.gradient(::MyDummyObjective2, x, p)
    return [x[1]^2 + x[2]^2 - 4, x[1] - x[2]]
end

function Cthonios.hessian(::MyDummyObjective2, x, p)
    return [
        2x[1] 2x[2];
        1. -1.
    ]
end

function test_newton_solver_func_2()
    objective = MyDummyObjective2()
    x = ones(2)
    p = nothing
    P = LUPreconditioner(objective, x, p)
    solver = Cthonios.NewtonSolver(objective, x, p)
    Cthonios.solve!(solver, objective, x, p, P)
    @test x[1] ≈ sqrt(2)
    @test x[2] ≈ sqrt(2)
end

function test_newton_solver()
    # test_newton_solver_func_1(false)
    # test_newton_solver_func_1(true)
    # test_newton_solver_func_1_bad_guess()
    test_newton_solver_func_2()
end

struct MyDummyObjective3 <: TestObjective
end

function Cthonios.value(::MyDummyObjective3, x, p)
    p = 1.
    return x[1] * (x[1] + 1) + 
           0.3 * x[2] * (x[2] - 0.2) + 
           0.2 * x[3] * (x[3] - 0.5) + 
           x[1] * x[1] * x[2] * x[2] + 
           p * x[1] * x[2] + 
           sin(x[1])

end

function Cthonios.gradient(o::MyDummyObjective3, x, p)
    return ForwardDiff.gradient(z -> Cthonios.value(o, z, p), x)
end

function Cthonios.hessian(o::MyDummyObjective3, x, p)
    return ForwardDiff.hessian(z -> Cthonios.value(o, z, p), x)
end

function Cthonios.hvp(o::MyDummyObjective3, x, v, p)
    return Cthonios.hessian(o, x, p) * v
end

function FiniteElementContainers.create_unknowns(::MyDummyObjective3)
    return zeros(3)
end

function test_trust_region_solver(verbose)
    objective = MyDummyObjective3()
    x = [2., 7., -1.]
    p = nothing
    preconditioner = NoPreconditioner(objective, x, p)
    solver = Cthonios.TrustRegionSolver(
        objective, x, p;
        verbose=verbose
    )
    Cthonios.solve!(solver, objective, x, p, preconditioner)
    g = Cthonios.gradient(objective, x, p)
    @test isapprox(norm(g), 0.; rtol=1.e7)
end

@testset "Newton solver" begin
    test_newton_solver()
end

@testset "Trust region solver" begin
    test_trust_region_solver(false)
    test_trust_region_solver(true)
end

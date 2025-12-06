function test_quasistatic_objective_constructors()
    obj = QuasiStaticObjective()
    @test obj.value == energy
    @test obj.gradient_u == residual
    @test obj.hessian_u == stiffness
end

function test_quasistatic_objective_cache_constructors(sim)
    objective = QuasiStaticObjective()
    cache, U, p = setup_caches(objective, sim)
end

function test_objectives()
    sim = sim_helper()
    test_quasistatic_objective_constructors()
    test_quasistatic_objective_cache_constructors(sim)
end

test_objectives()

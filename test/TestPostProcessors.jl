# uses some stuff from TestObjectives.jl, mainly sim_helper()
function test_post_processor_constructor_and_close(sim)
    objective = QuasiStaticObjective()
    cache, U, p = Cthonios.setup_caches(objective, sim)
    pp = Cthonios.PostProcessor(cache, U, p, sim)
    # do nothing for this test
    Cthonios.close(pp)
end

function test_post_processors()
    sim = sim_helper()
    test_post_processor_constructor_and_close(sim)
end

test_post_processors()

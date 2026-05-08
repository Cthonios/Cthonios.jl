function neohookean_props()
    props = Dict(
        "density"       => 1.0,
        "bulk modulus"  => 1000.,
        "shear modulus" => 1.0
    )
    return props
end

function test_sm_material_model_input_neohookean()
    props = neohookean_props()
    model = NeoHookean()
    formulation = PlaneStrain()
    physics = SolidMechanics(formulation, model)

    props = FiniteElementContainers.create_properties(physics, props)
    @test props[1] ≈ 1.0
    @test props[2] ≈ 1000.
    @test props[3] ≈ 1.

    state = FiniteElementContainers.create_initial_state(physics)
    @test state == SVector{0, Float64}()
end

test_sm_material_model_input_neohookean()

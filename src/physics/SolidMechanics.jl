struct SolidMechanics{NF, NP, NS, Form, Model} <: AbstractPhysics{NF, NP, NS}
    formulation::Form
    constitutive_model::Model
end

function SolidMechanics(formulation, model)
    NF = num_dimensions(formulation)
    NP = ConstitutiveModels.num_properties(model)
    NS = ConstitutiveModels.num_state_variables(model)
    return SolidMechanics{NF, NP, NS, typeof(formulation), typeof(model)}(
        formulation, model
    )
end

function FiniteElementContainers.create_properties(physics::SolidMechanics, inputs)
    density = inputs["density"]
    mat_model_props = ConstitutiveModels.initialize_props(physics.constitutive_model, inputs)
    mat_model_props = Array(mat_model_props)
    return pushfirst!(mat_model_props, density)
end

function FiniteElementContainers.create_initial_state(physics::SolidMechanics)
    return ConstitutiveModels.initialize_state(physics.constitutive_model)
end

function setup_function(fspace::FunctionSpace, ::SolidMechanics)
    return VectorFunction(fspace, :displ)
end

@inline function FiniteElementContainers.energy(
    physics::SolidMechanics, interps, x_el, t, dt, u_el, u_el_old, state_old_q, state_new_q, props_el

)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
  
    # kinematics
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    mat_props = @views props_el[2:end]

    # constitutive
    θ = 0.0 # TODO
    ψ_q = ConstitutiveModels.helmholtz_free_energy(
        physics.constitutive_model, mat_props, dt, ∇u_q, θ, state_old_q, state_new_q
    )  
    return JxW * ψ_q
end

# NOTE this does not integrate things
@inline function general_material_qoi(
    func,
    physics::SolidMechanics, interps, x_el, t, dt, u_el, u_el_old, state_old_q, state_new_q, props_el
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
  
    # kinematics
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    mat_props = @views props_el[2:end]

    # constitutive
    θ = 0.0 # TODO
    ψ_q = func(
        physics.constitutive_model, mat_props, dt, ∇u_q, θ, state_old_q, state_new_q
    )  
    return ψ_q
end

@inline function general_integrated_material_qoi(
    func,
    physics::SolidMechanics, interps, x_el, t, dt, u_el, u_el_old, state_old_q, state_new_q, props_el
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
  
    # kinematics
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    mat_props = @views props_el[2:end]

    # constitutive
    θ = 0.0 # TODO
    ψ_q = func(
        physics.constitutive_model, mat_props, dt, ∇u_q, θ, state_old_q, state_new_q
    )  
    return JxW * ψ_q
end

@inline function FiniteElementContainers.residual(
    physics::SolidMechanics, interps, x_el, t, dt, u_el, u_el_old, state_old_q, state_new_q, props_el
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
  
    # kinematics
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    mat_props = @views props_el[2:end]

    # constitutive
    θ = 0.0 # TODO
    P_q = ConstitutiveModels.pk1_stress(
        physics.constitutive_model, mat_props, dt, ∇u_q, θ, state_old_q, state_new_q
    )    
    # turn into voigt notation
    P_q = extract_stress(physics.formulation, P_q)
    G_q = discrete_gradient(physics.formulation, ∇N_X)
    f_q = G_q * P_q
    return JxW * f_q[:]
end

@inline function FiniteElementContainers.stiffness(  
    physics::SolidMechanics, interps, x_el, t, dt, u_el, u_el_old, state_old_q, state_new_q, props_el
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
    
    # kinematics
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    mat_props = @views props_el[2:end]

    # constitutive
    θ = 0. # TODO
    A_q = ConstitutiveModels.material_tangent(
        physics.constitutive_model, mat_props, dt, ∇u_q, θ, state_old_q, state_new_q
    )
    # turn into voigt notation
    K_q = extract_stiffness(physics.formulation, A_q)
    G_q = discrete_gradient(physics.formulation, ∇N_X)
    return JxW * G_q * K_q * G_q'
end

@inline function _kinetic_energy(
    physics::SolidMechanics, interps, x_el, t, dt, v_el, v_el_old, state_old_q, state_new_q, props_el

)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    v_q = interpolate_field_values(physics, interps, v_el)

    # TODO
    rho = props_el[1]
    return 0.5 * JxW * rho * dot(v_q, v_q)
end

@inline function kinetic_energy(
    physics::SolidMechanics, interps, x_el, t, dt, v_el, v_el_old, state_old_q, state_new_q, props_el
)
    return _kinetic_energy(physics, interps, v_el, x_el, state_old_q, props_el, t, dt)
end

function mass(
    physics::SolidMechanics, interps, x_el, t, dt, v_el, v_el_old, state_old_q, state_new_q, props_el
)
    return ForwardDiff.hessian(z -> _kinetic_energy(
        physics, interps, x_el, t, dt, z, v_el_old, state_old_q, state_new_q, props_el
    ), v_el)
end

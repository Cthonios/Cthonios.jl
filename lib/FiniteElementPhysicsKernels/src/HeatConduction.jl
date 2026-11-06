struct HeatConduction <: AbstractPhysics{1, 3, 0}
    constitutive_model::ConstitutiveModels.Thermal
end

function HeatConduction()
    return HeatConduction(ConstitutiveModels.Thermal(
        ConstitutiveModels.SimpleHeatCapacity(),
        FouriersLaw()
    ))
end

function FiniteElementContainers.create_properties(physics::HeatConduction, inputs)
    return ConstitutiveModels.initialize_props(physics.constitutive_model, inputs)
end

function setup_function(fspace::FunctionSpace, ::HeatConduction)
    return ScalarFunction(fspace, "temperature")
end

@inline function FiniteElementContainers.residual!(
    storage, e,
    physics::HeatConduction, t, dt, 
    props_el, state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    u_q, ∇u_q = interpolate_field_values_and_gradients(physics, interps, u_el)
    u_q_old = interpolate_field_values(physics, interps, u_el_old)
    dudt = (u_q - u_q_old) / dt
    form = GeneralFormulation{size(X_q, 1), num_fields(physics)}()
    ρ, θ_0, c, k = props_el
    scatter_with_values!(storage, form, e, conn, N, ρ * c * JxW * dudt)
    scatter_with_gradients!(storage, form, e, conn, ∇N_X, k * JxW * ∇u_q)
end

@inline function FiniteElementContainers.stiffness!(
    storage, e,
    physics::HeatConduction, t, dt,
    props_el, state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    form = GeneralFormulation{size(X_q, 1), num_fields(physics)}()
    ρ, θ_0, c, k = props_el
    scatter_with_values_and_values!(storage, form, e, conn, N, ρ * c * JxW)
    scatter_with_gradients_and_gradients!(storage, form, e, conn, ∇N_X, k * JxW)
end

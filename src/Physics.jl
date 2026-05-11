#############################################################################
# Material output helpers
#############################################################################
abstract type AbstractMaterialOutput end

struct StandardMaterialOutput{
    RT <: Number
} <: AbstractMaterialOutput
    algorithmic_energy::RT
    cauchy_stress::SymmetricTensor{2, 3, RT, 6}
    displacement_gradient::Tensor{2, 3, RT, 9}
    # state_variables::Svector{NS, RT}
end

@inline function standard_material_output(
    physics, interps, x_el, t, dt, u_el, u_el_old, state_old_q, state_new_q, props_el
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
    
    # kinematics
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    # constitutive
    θ = 0.0 # TODO
    ψ_q = ConstitutiveModels.helmholtz_free_energy(
        physics.constitutive_model, props_el, dt, ∇u_q, θ, state_old_q, state_new_q
    )  

    σ_q = ConstitutiveModels.cauchy_stress(
        physics.constitutive_model, props_el, dt, ∇u_q, θ, state_old_q, state_new_q
    )
    σ_q = σ_q |> symmetric
    return StandardMaterialOutput(ψ_q, σ_q, ∇u_q)
end

function update_material_output!(
    mat_output::L2Field,
    objective_cache, U, p
)
    update_field_unknowns!(p.field, assembler(objective_cache).dof, U)
    FiniteElementContainers.assemble_quadrature_quantity!(
        mat_output, nothing, assembler(objective_cache).dof,
        standard_material_output,
        U, p,
        FiniteElementContainers.AssembledStruct{StandardMaterialOutput}()
    )
    return nothing
end

function create_material_output(obj_cache, V_q, ::Type{T}) where T <: AbstractMaterialOutput
    fspace = FiniteElementContainers.function_space(assembler(obj_cache).dof)

    mat_outputs = StructArray{T}[]
    mat_vars = FiniteElementContainers.AbstractFunction[]

    mat_vars = (
        ScalarFunction(V_q, "algorithmic_energy"),
        SymmetricTensorFunction(V_q, "cauchy_stress"),
        TensorFunction(V_q, "displacement_gradient")
    )

    mat_outputs = L2Field(undef, T, 1, FiniteElementContainers.block_quadrature_sizes(fspace))

    return mat_outputs, mat_vars
end

#############################################################################
# Solid mechanics methods
#############################################################################
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
    return VectorFunction(fspace, "displ")
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
    # P_q = tovoigt(SVector, P_q)
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
    # K_q = tovoigt(SMatrix, A_q)
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

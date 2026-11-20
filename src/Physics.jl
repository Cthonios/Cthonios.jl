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
        physics.constitutive_model, props_el, state_old_q, state_new_q, dt, ∇u_q, θ
    )  

    σ_q = ConstitutiveModels.cauchy_stress(
        physics.constitutive_model, props_el, state_old_q, state_new_q, dt, ∇u_q, θ
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

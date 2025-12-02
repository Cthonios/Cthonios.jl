abstract type AbstractMaterialOutput end

struct StandardMaterialOutput{
    RT <: Number
} <: AbstractMaterialOutput
    algorithmic_energy::RT
    cauchy_stress::SymmetricTensor{2, 3, RT, 6}
    displacement_gradient::Tensor{2, 3, RT, 9}
end

@inline function standard_material_output(
    physics, interps, x_el, t, dt, u_el, u_el_old, state_old_q, state_new_q, props_el
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
    
    # kinematics
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)
    # F_q = ∇u_q + one(∇u_q)

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
    mat_output::NamedTuple,
    objective_cache
)
    T = eltype(values(mat_output)[1])
    U = objective_cache.solution
    p = objective_cache.parameters
    FiniteElementContainers.assemble_quadrature_quantity!(
        mat_output, assembler(objective_cache).dof,
        standard_material_output,
        U, p
    )
    return nothing
end

function create_material_output(obj_cache, V_q, ::Type{T}) where T <: AbstractMaterialOutput
    fspace = FiniteElementContainers.function_space(assembler(obj_cache).dof)
    syms = keys(fspace.elem_conns)

    mat_outputs = StructArray{T}[]
    mat_vars = FiniteElementContainers.AbstractFunction[]

    mat_vars = (
        ScalarFunction(V_q, :algorithmic_energy),
        SymmetricTensorFunction(V_q, :cauchy_stress),
        TensorFunction(V_q, :displacement_gradient)
    )

    for (n, sym) in enumerate(syms)
        conn = getproperty(fspace.elem_conns, sym)
        ref_fe = getproperty(fspace.ref_fes, sym)
        nq = num_quadrature_points(ref_fe)
        ne = size(conn, 2)
        temp_output = StructArray{T}(undef, nq, ne)
        push!(mat_outputs, temp_output)
    end

    mat_outputs = NamedTuple{syms}(tuple(mat_outputs...))

    # V = FunctionSpace()
    return mat_outputs, mat_vars
end

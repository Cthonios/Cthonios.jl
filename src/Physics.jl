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
# Heat conduction methods
############################################################################
struct HeatConduction <: AbstractPhysics{1, 3, 0}
    constitutive_model::ConstitutiveModels.FouriersLaw
end

function HeatConduction()
    return HeatConduction(ConstitutiveModels.FouriersLaw())
end

function FiniteElementContainers.create_properties(physics::HeatConduction, inputs)
    density = inputs["density"]::Float64
    cp = inputs["heat capacity"]
    k = inputs["thermal conductivity"]
    # mat_model_props = ConstitutiveModels.initialize_props(physics.constitutive_model, inputs)
    # mat_model_props = Array(mat_model_props)
    # return mat_model_props
    return [density, cp, k]
end

function setup_function(fspace::FunctionSpace, ::HeatConduction)
    return ScalarFunction(fspace, "temperature")
end

@inline function FiniteElementContainers.residual!(
    storage, e,
    physics::HeatConduction, t, dt, props_el, 
    state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    u_q, ∇u_q = interpolate_field_values_and_gradients(physics, interps, u_el)
    u_q_old = interpolate_field_values(physics, interps, u_el_old)
    dudt = (u_q - u_q_old) / dt
    # R_q = ∇u_q * ∇N_X'
    form = GeneralFormulation{size(X_q, 1), num_fields(physics)}()
    density, cp, k = props_el
    scatter_with_values!(storage, form, e, conn, N, density * cp * JxW * dudt)
    scatter_with_gradients!(storage, form, e, conn, ∇N_X, k * JxW * ∇u_q)
end

@inline function FiniteElementContainers.stiffness!(
    storage, e,
    physics::HeatConduction, t, dt, props_el, 
    state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    # ∇u_q = interpolate_field_gradients(physics, interps, u_el)
    # R_q = ∇u_q * ∇N_X' - N' * physics.func(X_q, 0.0)
    form = GeneralFormulation{size(X_q, 1), num_fields(physics)}()
    density, cp, k = props_el
    scatter_with_values_and_values!(storage, form, e, conn, N, density * cp * JxW)
    scatter_with_gradients_and_gradients!(storage, form, e, conn, ∇N_X, k * JxW)
end

#############################################################################
# Solid mechanics methods
############################################################################
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
    density = inputs["density"]::Float64
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

function characteristic_element_length(
    physics::SolidMechanics, interps, x_el,
    t, dt, u_el, u_el_old,
    state_old_q, state_new_q, props_el,
)
    x_cur = x_el + u_el
    ndim = 3
    nnpe = length(x_cur) ÷ ndim
    T = eltype(x_cur)
    cx = cy = cz = zero(T)
    for i in 1:nnpe
        cx += x_cur[(i-1)*ndim + 1]
        cy += x_cur[(i-1)*ndim + 2]
        cz += x_cur[(i-1)*ndim + 3]
    end
    cx /= nnpe; cy /= nnpe; cz /= nnpe
    total = zero(T)
    for i in 1:nnpe
        dx = x_cur[(i-1)*ndim + 1] - cx
        dy = x_cur[(i-1)*ndim + 2] - cy
        dz = x_cur[(i-1)*ndim + 3] - cz
        total += sqrt(dx*dx + dy*dy + dz*dz)
    end
    return 2 * total / nnpe
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

@inline function FiniteElementContainers.mass(
    physics::SolidMechanics,
    interps, x_el,
    t, dt,
    u_el, u_el_old,
    state_old_q, state_new_q,
    props_el,
)
    cell = map_interpolants(interps, x_el)
    (; N, JxW) = cell

    # Build element mass matrix in interleaved DOF ordering:
    #   M_el[3*(n-1)+d, 3*(m-1)+d'] = δ(d,d') * N[n] * N[m]
    # i.e. kron(N*N', I_3).  The FEC assembly infrastructure expects
    # rows/cols in the same interleaved order as discrete_gradient, so
    # "N_vec * N_vec'" with a block-ordered N_vec would be wrong.
    N_nodes = size(N, 1)
    NDOF    = 3 * N_nodes
    tup = zeros(SVector{NDOF * NDOF, eltype(N)})
    for n in 1:N_nodes
        for m in 1:N_nodes
            Nnm = N[n] * N[m]
            for d in 1:3
                r = 3 * (n - 1) + d
                c = 3 * (m - 1) + d
                linear_idx = r + NDOF * (c - 1)   # column-major flat index
                tup = setindex(tup, Nnm, linear_idx)
            end
        end
    end
    ρ = props_el[1]
    M_el = SMatrix{NDOF, NDOF, eltype(N), NDOF * NDOF}(tup.data)
    return JxW * ρ * M_el
end

@inline function FiniteElementContainers.mass!(
    storage, e,
    physics::SolidMechanics, t, dt, props_el, 
    state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old
)
    interps = map_interpolants(interps, x_el)
    (; X_q, N, ∇N_X, JxW) = interps
    # ∇u_q = interpolate_field_gradients(physics, interps, u_el)
  
    # kinematics
    # ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    # mat_props = @views props_el[2:end]
    ρ = props_el[1]

    # constitutive
    # θ = 0. # TODO
    # A_q = ConstitutiveModels.material_tangent(
    #     physics.constitutive_model, mat_props, dt, ∇u_q, θ, state_old_q, state_new_q
    # )
    scatter_with_values_and_values!(
      storage, physics.formulation, e, conn, N, JxW * ρ
    )
    return nothing
end

@inline function FiniteElementContainers.mass_action(
    physics::SolidMechanics,
    interps, x_el,
    t, dt,
    u_el, u_el_old, v_el,
    state_old_q, state_new_q,
    props_el,
)
    cell = map_interpolants(interps, x_el)
    (; N, JxW) = cell
    ρ = props_el[1]
    # Correct M·v in interleaved DOF ordering:
    #   (M·v)[3*(n-1)+d] = N[n] * Σ_m N[m] * v_el[3*(m-1)+d]
    # i.e. per-direction dot products, NOT a single dot over all DOFs.
    N_nodes = size(N, 1)
    s1 = sum(N[m] * v_el[3 * (m - 1) + 1] for m in 1:N_nodes)
    s2 = sum(N[m] * v_el[3 * (m - 1) + 2] for m in 1:N_nodes)
    s3 = sum(N[m] * v_el[3 * (m - 1) + 3] for m in 1:N_nodes)
    tup = zeros(SVector{3 * N_nodes, eltype(v_el)})
    for n in 1:N_nodes
        k = 3 * (n - 1)
        tup = setindex(tup, N[n] * s1, k + 1)
        tup = setindex(tup, N[n] * s2, k + 2)
        tup = setindex(tup, N[n] * s3, k + 3)
    end
    return JxW * ρ * tup
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

@inline function FiniteElementContainers.residual!(
    storage, e,
    physics::SolidMechanics, t, dt, props_el, 
    state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old
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
    scatter_with_gradients!(storage, physics.formulation, e, conn, ∇N_X, JxW * P_q)
    return nothing
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

@inline function FiniteElementContainers.stiffness!(
    storage, e,
    physics::SolidMechanics, t, dt, props_el, 
    state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old
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
    scatter_with_gradients_and_gradients!(
      storage, physics.formulation, e, conn, ∇N_X, JxW * A_q
    )
    return nothing
end

@inline function FiniteElementContainers.stiffness_action(
    physics::SolidMechanics, interps, x_el, t, dt, u_el, u_el_old, v_el, state_old_q, state_new_q, props_el
)
    interps = map_interpolants(interps, x_el)
    (; ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
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
    return JxW * G_q * (K_q * (G_q' * v_el))
end

@inline function FiniteElementContainers.stiffness_action!(
    storage, e,
    physics::SolidMechanics, t, dt, props_el, 
    state_old_q, state_new_q,
    conn, interps, x_el, u_el, u_el_old, v_el
  )
    interps = map_interpolants(interps, x_el)
    (; ∇N_X, JxW) = interps
    ∇u_q = interpolate_field_gradients(physics, interps, u_el)
    ∇u_q = modify_field_gradients(physics.formulation, ∇u_q)

    mat_props = @views props_el[2:end]

    # constitutive
    θ = 0. # TODO
    A_q = ConstitutiveModels.material_tangent(
        physics.constitutive_model, mat_props, dt, ∇u_q, θ, state_old_q, state_new_q
    )
    scatter_with_gradients_and_gradients!(storage, physics.formulation, e, conn, ∇N_X, JxW * A_q, v_el)
    return nothing
end

struct TotalLagrangeSection{
  FS, Form, Mod, Props
} <: Section{FS, Form, Mod}
  fspace::FS
  formulation::Form
  model::Mod
  props::Props
end

function TotalLagrangeSection(
  fspace::F, formulation::Form,
  model::M, props::P, initial_state::S
) where {
  F <: FunctionSpace, Form <: AbstractMechanicsFormulation,
  M <: MechanicalModel, P <: AbstractVector{<:Number}, S
}

  return TotalLagrangeSection{
    typeof(fspace), typeof(formulation),
    typeof(model), typeof(props)
  }(fspace, formulation, model, props)
end

function Base.show(io::IO, section::TotalLagrangeSection)
  print(io, "          TotalLagrangeSection\n",
        "            Function space            = $(typeof(section.fspace).name.name)\n",
        "            Formulation               = $(section.formulation)\n",
        "            Material model            = $(section.model)\n",
        "            Material model properties = $(section.props)\n")
end

"""
Quadrature point energy method
"""
function energy(section::Section, U_el, X_el, q)
  # unpack some stuff, TODO state vars need to not be Hardcoded
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  # ∇N_X, JxW, _ = interpolants(section.fspace.ref_fe, section.formulation, X_el, q)

  ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
  J    = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'
  JxW  = det(J) * ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe, q)

  ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)
  ψ_q = ConstitutiveModels.energy(model, props, F_q, state)
  return JxW * ψ_q
end

"""
Quadrature point residual method
"""
function residual(section::Section, U_el, X_el, q)
  # unpack some stuff, TODO state vars need to not be Hardcoded
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  # ∇N_X, JxW, _ = interpolants(section.fspace.ref_fe, section.formulation, X_el, q)

  ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
  J    = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'
  JxW  = det(J) * ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe, q)

  ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)
  P_q = ConstitutiveModels.pk1_stress(model, props, F_q, state)
  P_v = FiniteElementContainers.extract_stress(section.formulation, P_q) 
  G = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)
  return JxW * G * P_v
end

"""
Quadrature point hvp method
"""
function hvp(section::Section, U_el, V_el, X_el, q)
  # model, props = section.model, section.props
  # state        = SVector{0, Float64}()
  # ∇N_X, JxW, _ = interpolants(section.fspace.ref_fe, section.formulation, X_el, q)

  # ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
  # J    = X_el * ∇N_ξ
  # J_inv = inv(J)
  # ∇N_X = (J_inv * ∇N_ξ')'
  # JxW  = det(J) * ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe, q)

  # ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
  # F_q = ∇u_q + one(∇u_q)

  # # A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)
  # # A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)
  # # G    = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)

  # P_q = ConstitutiveModels.pk1_stress(model, props, F_q, state)
  # P_v = FiniteElementContainers.extract_stress(section.formulation, P_q) 
  # G = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)
  # P_new = G * P_v
  # # @show size(P_v) size(G) size(V_el)

  # grad = ForwardDiff.gradient(x -> dot(x, V_el), P_new)
  # display(grad)
  # # @assert false
  # A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)
  # A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)
  # display(G * A_m * G' * V_el)

  # @show energy(section, U_el, X_el, q)
  # R = ForwardDiff.gradient(x -> energy(section, x, X_el, q), U_el)
  # @show R

  Hv = ForwardDiff.gradient(x -> dot(ForwardDiff.gradient(y -> energy(section, y, X_el, q), x), V_el), U_el)
  # @show Hv
  # @assert false
  # return JxW * dot(ForwardDiff.gradient(x -> dot(x, V_el), P_new)[1], G)
  # return JxW * G * A_m * G' * V_el
  return vec(Hv)
end

"""
Quadrature point stiffness method
"""
function stiffness(section::Section, U_el, X_el, q)
  # unpack some stuff, TODO state vars need to not be Hardcoded
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  # ∇N_X, JxW, _ = interpolants(section.fspace.ref_fe, section.formulation, X_el, q)

  ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
  J    = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X = (J_inv * ∇N_ξ')'
  JxW  = det(J) * ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe, q)

  ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
  F_q = ∇u_q + one(∇u_q)
  A_q = ConstitutiveModels
  A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)
  A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)
  G    = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)
  return JxW * G * A_m * G'
end

"""
Top level energy method for a section
Takes in nodal fields for the fields and coordinates
"""
function energy(section::Section, U::V1, X::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}

  # ND, NN = num_dimensions(section), num_nodes_per_element(section)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  W      = 0.0
  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)
    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    for q in 1:FiniteElementContainers.num_q_points(section)
      W = W + energy(section, U_el, X_el, q)
    end
  end
  W
end

"""
Top level residual method for a section
"""
function residual!(asm::Assembler, section::Section, U::V1, X::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)
    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    R_el = zero(SVector{ND * NN, eltype(U)})
    for q in 1:FiniteElementContainers.num_q_points(section)
      R_el = R_el + residual(section, U_el, X_el, q)
    end
    assemble!(asm, R_el, conn)
  end
end

"""
Top level hvp method for a section
"""
function hvp!(asm::Assembler, section::Section, U::V1, V::V2, X::V3) where {V1 <: AbstractArray, V2 <: AbstractArray, V3 <: AbstractArray}
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)
    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    # V_el = SVector{ND * NN, eltype(X)}(@views V[conn])
    V_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views V[conn])
    Av_el = zero(SVector{ND * NN, eltype(U)})
    for q in 1:FiniteElementContainers.num_q_points(section)
      Hv = SVector{ND * NN, eltype(Av_el)}(hvp(section, U_el, V_el, X_el, q) |> vec)
      # Av_el = Av_el + hvp(section, U_el, V_el, X_el, q)
      Av_el = Av_el + Hv
    end
    assemble!(asm, Av_el, conn)
  end
end

"""
Top level stiffness method for a section
"""
function stiffness!(asm::Assembler, section::Section, U::V1, X::V2, sec_id::Int) where {V1 <: AbstractArray, V2 <: AbstractArray}
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  for e in 1:num_elements(section)
    conn = dof_connectivity(section, e)
    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views U[conn])
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views X[conn])
    K_el = zero(SMatrix{ND * NN, ND * NN, Float64, ND * NN * ND * NN})
    for q in 1:FiniteElementContainers.num_q_points(section)
      K_el = K_el + stiffness(section, U_el, X_el, q)
    end
    assemble!(asm, K_el, sec_id, e)
  end
end

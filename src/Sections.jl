abstract type AbstractSection <: AbstractCthoniosType end

struct Section{
  FSpace   <: FunctionSpace,
  Form     <: AbstractMechanicsFormulation,
  Centroid <: ReferenceFE,
  Model    <: ConstitutiveModel,
  Props    <: AbstractArray{<:Number, 1}
  # State  <: QuadratureField
}
  fspace::FSpace
  formulation::Form
  centroid_ref_fe::Centroid
  model::Model
  props::Props
  # state_old::State
  # state_new::State
end

num_dimensions(sec::Section)           = FiniteElementContainers.num_dimensions(sec.fspace)
num_elements(sec::Section)             = FiniteElementContainers.num_elements(sec.fspace)
num_nodes_per_element(sec::Section)    = FiniteElementContainers.num_nodes_per_element(sec.fspace)
num_q_points(sec::Section)             = FiniteElementContainers.num_q_points(sec.fspace)

connectivity(sec::Section)             = FiniteElementContainers.connectivity(sec.fspace)
connectivity(sec::Section, e::Int)     = FiniteElementContainers.connectivity(sec.fspace, e)
dof_connectivity(sec::Section)         = FiniteElementContainers.dof_connectivity(sec.fspace)
dof_connectivity(sec::Section, e::Int) = FiniteElementContainers.dof_connectivity(sec.fspace, e)

element_level_fields(sec::Section, U) = 
FiniteElementContainers.element_level_fields(sec.fspace, U)

element_level_fields(sec::Section, U, e::Int) =
FiniteElementContainers.element_level_fields(sec.fspace, U, e)

# shape_function_values(sec::Section, X_el, q::Int, ::Int) = 
# ReferenceFiniteElements.shape_function_values(sec.fspace.ref_fe, q)

# function shape_function_gradients(sec::Section, X_el)

# end

function interpolants(ref_fe::ReferenceFE, formulation, X_el, q::Int)
  ∇N_ξ  = ReferenceFiniteElements.shape_function_gradients(ref_fe, q)
  J     = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X  = (J_inv * ∇N_ξ')'
  JxW   = det(J) * ReferenceFiniteElements.quadrature_weight(ref_fe, q)
  G     = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return ∇N_X, JxW, G
end

# function modify_deformation_gradient(::Type{<:FiniteElementContainers.AbstractMechanicsFormulation}, F_q, F_c)
#   if section.formulation <: FiniteElementContainers.PlaneStrain
#     F_q = ((det(F_c) / det(F_q))^(1. / 2.)) * F_q
#   else
#     @assert false "unsupported formulation"
#   end
#   return F_q
# end

function Fbar_deformation_gradient(::FiniteElementContainers.PlaneStrain, F_q, F_c)
  F_q = ((det(F_c) / det(F_q))^(1. / 2.)) * F_q
  return F_q
end

function Fbar_stifness(::FiniteElementContainers.PlaneStrain, G, G_c, P_v, A_m)
  T   = eltype(P_v)
  TxI = SMatrix{4, 4, T, 16}((
    P_v[1],  P_v[2],  P_v[3],  P_v[4],
    zero(T), zero(T), zero(T), zero(T),
    zero(T), zero(T), zero(T), zero(T),
    P_v[1],  P_v[2],  P_v[3],  P_v[4],
  ))
  T   = eltype(A_m)
  A_II = SMatrix{4, 4, T, 16}((
    A_m[1, 1] + A_m[1, 4], A_m[2, 1] + A_m[2, 4], A_m[3, 1] + A_m[3, 4], A_m[4, 1] + A_m[4, 4],
    zero(T),               zero(T),               zero(T),               zero(T),
    zero(T),               zero(T),               zero(T),               zero(T),
    A_m[1, 1] + A_m[1, 4], A_m[2, 1] + A_m[2, 4], A_m[3, 1] + A_m[3, 4], A_m[4, 1] + A_m[4, 4]
  ))
  Q = 0.5 * A_II - 0.5 * TxI

  K_bar = G * Q * (G_c - G)'
  return K_bar
end

###################################################################

function Section(
  fspace::F, formulation::Form,
  model::M, props::P, initial_state::S; Fbar::Bool=false
) where {F    <: FunctionSpace, 
         Form <: AbstractMechanicsFormulation,
         M    <: ConstitutiveModels.MechanicalModel,
         P    <: AbstractArray{<:Number, 1},
         S    <: AbstractArray}

  # unpack some compile time constants
  # Q  = FiniteElementContainers.num_q_points(fspace)
  # E  = FiniteElementContainers.num_elements(fspace)
  # NF = length(initial_state)
  # NF        = 1

  # set up state from a simple initial state
  # state_old = QuadratureField{length(initial_state), Q, E, Vector, typeof(initial_state)}(undef)
  # state_new = QuadratureField{length(initial_state), Q, E, Vector, typeof(initial_state)}(undef)
  # state_old .= initial_state 
  # state_new .= initial_state

  ref_fe_type = ReferenceFE(@eval $(typeof(fspace.ref_fe.ref_fe_type).name.name)(Val(1)))


  return Section{
    typeof(fspace), typeof(formulation), typeof(ref_fe_type),
    typeof(model), typeof(props)#, typeof(state_old)  
  }(fspace, formulation, ref_fe_type, model, props)#, state_old, state_new)
end


I = one(Tensor{2, 3, Float64, 9})

# ref_fe_plane_strain_quad4 = ReferenceFE(Quad4(Val(1)))

# function F_bar(::Type{PlaneStrain}, F)

# end

# currently no state vars setup TODO
# function energy(section::Section, X::NodalField, U::NodalField)

#   ND, NN       = num_dimensions(section), num_nodes_per_element(section)
#   model, props = section.model, section.props
#   state        = SVector{0, Float64}()
#   W = 0.0
#   for e in 1:num_elements(section)
#     conn = connectivity(section, e)
#     X_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(X[:, conn]))
#     U_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(U[:, conn]))
#     for q in 1:num_q_points(section)
#       ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
#       J    = X_el * ∇N_ξ
#       J_inv = inv(J)
#       ∇N_X = (J_inv * ∇N_ξ')'
#       JxW  = det(J) * ReferenceFiniteElements.quadrature_weight(section.fspace.ref_fe, q)
#       ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
#       F_q  = ∇u_q + one(∇u_q)
#       W_q  = ConstitutiveModels.energy(model, props, F_q, state)
#       W    = W + JxW * W_q
#     end
#   end
#   W
# end

function residual!(R, section::Section, X::NodalField, U::NodalField)
  ND, NN       = num_dimensions(section), num_nodes_per_element(section)
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  # FiniteElementContainers.update_fields!(domain.U, domain.dof, Uu)
  for e in 1:num_elements(section)
    conn = connectivity(section, e)
    dof_conn = dof_connectivity(section, e)
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views vec(U[:, conn]))
    # R_el = zero(SMatrix{ND, NN, Float64, ND * NN})
    R_el = zero(SVector{ND * NN, eltype(U)}) # Unitful issue here
    for q in 1:num_q_points(section)
      ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
      J    = X_el * ∇N_ξ
      J_inv = inv(J)
      ∇N_X = (J_inv * ∇N_ξ')'
      JxW  = det(J) * ReferenceFiniteElements.quadrature_weight(section.fspace.ref_fe, q)
      ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
      F_q  = ∇u_q + one(∇u_q)
      P_q  = ConstitutiveModels.pk1_stress(model, props, F_q, state)
      P_v  = FiniteElementContainers.extract_stress(section.formulation, P_q) 
      G    = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)
      R_el   = R_el + JxW * G * P_v
      # R_el = R_el + JxW * P_q[1:2, 1:2] * ∇N_X' # Hardcoded to 2D problems right now
    end
    FiniteElementContainers.assemble_residual!(R, R_el, dof_conn)
  end
end

function stiffness!(K, section::Section, X::NodalField, U::NodalField)
  ND, NN       = num_dimensions(section), num_nodes_per_element(section)
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  for e in 1:num_elements(section)
    conn = connectivity(section, e)
    dof_conn = dof_connectivity(section, e)
    X_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(U[:, conn]))
    K_el = zero(SMatrix{ND * NN, ND * NN, Float64, ND * NN * ND * NN})
    for q in 1:num_q_points(section)
      ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
      J    = X_el * ∇N_ξ
      J_inv = inv(J)
      ∇N_X = (J_inv * ∇N_ξ')'
      JxW  = det(J) * ReferenceFiniteElements.quadrature_weight(section.fspace.ref_fe, q)
      ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
      F_q  = ∇u_q + one(∇u_q)
      A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)
      A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)
      G    = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)
      K_el = K_el + JxW * G * A_m * G'
    end
    FiniteElementContainers.assemble!(K, K_el, dof_conn)
  end
end

#####################################################################

function energy(section::Section, X_el, U_el, q)
  # unpack some stuff, TODO state vars need to not be Hardcoded
  model, props = section.model, section.props
  state        = SVector{0, Float64}()

  # interpolants, also get centroid element
  ∇N_X, JxW, _   = interpolants(section.fspace.ref_fe, section.formulation, X_el, q)
  # ∇N_X_c, _, G_c = interpolants(section.centroid_ref_fe, section.formulation, X_el, 1)

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
  # ∇u_c = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X_c)
  F_q  = ∇u_q + one(∇u_q)
  # F_c  = ∇u_c + one(∇u_c)

  # TODO wrap in if statement with Fbar as option
  # F_bar = F_q
  # F_bar = Fbar_deformation_gradient(section.formulation, F_q, F_c)
  F_bar = F_q

  # constitutive model
  ψ_q = ConstitutiveModels.energy(model, props, F_bar, state)

  return JxW * ψ_q
end

function residual_and_stiffness(section::Section, X_el, U_el, q)
  
  # unpack some stuff, TODO state vars need to not be Hardcoded
  model, props = section.model, section.props
  state        = SVector{0, Float64}()

  # interpolants, also get centroid element
  ∇N_X, JxW, G   = interpolants(section.fspace.ref_fe, section.formulation, X_el, q)
  ∇N_X_c, _, G_c = interpolants(section.centroid_ref_fe, section.formulation, X_el, 1)

  # kinematics
  ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
  ∇u_c = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X_c)
  F_q  = ∇u_q + one(∇u_q)
  F_c  = ∇u_c + one(∇u_c)

  # TODO wrap in if statement with Fbar as option
  # F_bar = F_q
  # F_bar = Fbar_deformation_gradient(section.formulation, F_q, F_c)
  F_bar = F_q

  # constitutive model
  P_q  = ConstitutiveModels.pk1_stress(model, props, F_bar, state)
  A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)

  # vectorization
  P_v  = FiniteElementContainers.extract_stress(section.formulation, P_q)
  A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)

  # accumulation to element level
  R_q = JxW * G * P_v
  # K_q = JxW * G * A_m * G'

  K     = G * A_m * G'

  # TODO add if statement here for Fbar as an option
  K_bar = Fbar_stifness(section.formulation, G, G_c, P_v, A_m)

  # K_q   = JxW * (K + K_bar)
  K_q = JxW * K

  return R_q, K_q
end

# function energy(section::Section, X::NodalField, U::NodalField)
function energy(section::Section, X::V1, U::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}

  ND, NN = num_dimensions(section), num_nodes_per_element(section)
  W      = 0.0
  for e in 1:num_elements(section)
    conn = connectivity(section, e)
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views vec(U[:, conn]))
    for q in 1:num_q_points(section)
      W = W + energy(section, X_el, U_el, q)
    end
  end
  W
end

# I don't think you need to touch this one anymore
function assemble!(assembler::StaticAssembler, section::Section, X::NodalField, U::NodalField)
  ND, NN       = num_dimensions(section), num_nodes_per_element(section)
  for e in 1:num_elements(section)

    # TODO eliminate one of below somehow
    conn = connectivity(section, e)
    dof_conn = dof_connectivity(section, e)
    X_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(U[:, conn]))

    R_el = zero(SVector{ND * NN, Float64})
    K_el = zero(SMatrix{ND * NN, ND * NN, Float64, ND * NN * ND * NN})

    # quadrature point loop
    for q in 1:num_q_points(section)
      # TODO put hook here for Fbar
      R_q, K_q = residual_and_stiffness(section, X_el, U_el, q)

      R_el += R_q
      K_el += K_q
    end
    FiniteElementContainers.assemble!(assembler.R, R_el, dof_conn)
    FiniteElementContainers.assemble!(assembler.K, K_el, dof_conn)
  end
end

##################################################

function jvp!(y, section::Section, X::NodalField, U::NodalField, v::NodalField)
  ND, NN       = num_dimensions(section), num_nodes_per_element(section)
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  for e in 1:num_elements(section)

    # TODO eliminate one of below somehow
    conn = connectivity(section, e)
    dof_conn = dof_connectivity(section, e)
    X_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(U[:, conn]))
    # v_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(v[:, conn]))
    # v_el = SVector{ND * NN, Float64}(@views vec(v[:, conn]))
    v_el = SVector{ND * NN, Float64}(@views v[dof_conn])

    # R_el = zero(SVector{ND * NN, Float64})
    # K_el = zero(SMatrix{ND * NN, ND * NN, Float64, ND * NN * ND * NN})
    y_el = zero(SVector{ND * NN, Float64})

    # quadrature point loop
    for q in 1:num_q_points(section)
      # interpolants
      ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
      J    = X_el * ∇N_ξ
      J_inv = inv(J)
      ∇N_X = (J_inv * ∇N_ξ')'
      JxW  = det(J) * ReferenceFiniteElements.quadrature_weight(section.fspace.ref_fe, q)
      G    = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)

      # kinematics
      ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
      F_q  = ∇u_q + one(∇u_q)

      # constitutive model
      # P_q  = Tensors.gradient(z -> ConstitutiveModels.energy(model, props, z, state), F_q)
      # P_q  = ConstitutiveModels.pk1_stress(model, props, F_q, state)
      A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)
      # P_v  = FiniteElementContainers.extract_stress(section.formulation, P_q)
      A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)
      # R_el = R_el + JxW * G * P_v
      # K_el = K_el + JxW * G * A_m * G'
      y_el = y_el + JxW * G * A_m * G' * v_el
    end
    # FiniteElementContainers.assemble!(assembler.R, R_el, dof_conn)
    # FiniteElementContainers.assemble!(assembler.K, K_el, dof_conn)
    FiniteElementContainers.assemble!(y, y_el, dof_conn)
  end
end
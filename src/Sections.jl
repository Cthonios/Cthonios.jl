abstract type AbstractSection <: AbstractCthoniosType end

struct Section{
  FSpace <: FunctionSpace,
  Form   <: AbstractMechanicsFormulation,
  Model  <: ConstitutiveModel,
  Props  <: AbstractArray{<:Number, 1}
  # State  <: QuadratureField
}
  fspace::FSpace
  formulation::Form
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

function Section(
  fspace::F, formulation::Form,
  model::M, props::P, initial_state::S
) where {F    <: FunctionSpace, 
         Form <: AbstractMechanicsFormulation,
         M    <: ConstitutiveModels.MechanicalModel,
         P    <: AbstractArray{<:Number, 1},
         S    <: AbstractArray}

  # unpack some compile time constants
  Q  = FiniteElementContainers.num_q_points(fspace)
  E  = FiniteElementContainers.num_elements(fspace)
  NF = length(initial_state)
  # NF        = 1

  @show typeof(initial_state)
  # set up state from a simple initial state
  # state_old = QuadratureField{length(initial_state), Q, E, Vector, typeof(initial_state)}(undef)
  # state_new = QuadratureField{length(initial_state), Q, E, Vector, typeof(initial_state)}(undef)
  # state_old .= initial_state 
  # state_new .= initial_state

  return Section{
    typeof(fspace), typeof(formulation), typeof(model),
    typeof(props)#, typeof(state_old)  
  }(fspace, formulation, model, props)#, state_old, state_new)
end


I = one(Tensor{2, 3, Float64, 9})

# currently no state vars setup TODO
function energy(section::Section, X::NodalField, U::NodalField)

  ND, NN       = num_dimensions(section), num_nodes_per_element(section)
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  W = 0.0
  for e in 1:num_elements(section)
    conn = connectivity(section, e)
    X_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(U[:, conn]))
    for q in 1:num_q_points(section)
      ∇N_ξ = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
      J    = X_el * ∇N_ξ
      J_inv = inv(J)
      ∇N_X = (J_inv * ∇N_ξ')'
      JxW  = det(J) * ReferenceFiniteElements.quadrature_weight(section.fspace.ref_fe, q)
      ∇u_q = FiniteElementContainers.modify_field_gradients(section.formulation, U_el * ∇N_X)
      F_q  = ∇u_q + one(∇u_q)
      W_q  = ConstitutiveModels.energy(model, props, F_q, state)
      W    = W + JxW * W_q
    end
  end
  W
end

function residual!(R, section::Section, X::NodalField, U::NodalField)
  ND, NN       = num_dimensions(section), num_nodes_per_element(section)
  model, props = section.model, section.props
  state        = SVector{0, Float64}()
  for e in 1:num_elements(section)
    conn = connectivity(section, e)
    dof_conn = dof_connectivity(section, e)
    X_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, Float64, ND * NN}(@views vec(U[:, conn]))
    # R_el = zero(SMatrix{ND, NN, Float64, ND * NN})
    R_el = zero(SVector{ND * NN, Float64})
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
    assemble!(R, R_el, dof_conn)
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
      # A_q  = ConstitutiveModels.material_tangent(model, props, F_q, state)
      A_q  = Tensors.hessian(z -> ConstitutiveModels.energy(model, props, z, state), F_q)
      A_m  = FiniteElementContainers.extract_stiffness(section.formulation, A_q)
      G    = FiniteElementContainers.discrete_gradient(section.formulation, ∇N_X)
      # K_el = K_el + JxW * G * A_m
      K_el = K_el + JxW * G * A_m * G'
      # @show size(G)
    end
    assemble!(K, K_el, dof_conn)
  end
end

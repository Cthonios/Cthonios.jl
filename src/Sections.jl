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

function Base.show(io::IO, section::Section)
  print(io, "          Section\n",
        "            Function space            = $(section.fspace)\n",
        "            Formulation               = $(section.formulation)\n",
        "            Material model            = $(section.model)\n",
        "            Material model properties = $(section.model)\n")
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

function interpolants(ref_fe::ReferenceFE, formulation, X_el, q::Int)
  ∇N_ξ  = ReferenceFiniteElements.shape_function_gradients(ref_fe, q)
  J     = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X  = (J_inv * ∇N_ξ')'
  JxW   = det(J) * ReferenceFiniteElements.quadrature_weight(ref_fe, q)
  G     = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return ∇N_X, JxW, G
end

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
      # R_el = R_el + JxW * vec(P_q[1:2, 1:2] * ∇N_X') # Hardcoded to 2D problems right now
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
    X_el = SMatrix{ND, NN, eltype(X), ND * NN}(@views vec(X[:, conn]))
    U_el = SMatrix{ND, NN, eltype(U), ND * NN}(@views vec(U[:, conn]))
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

# Parsing
function read_materials(input_settings::D) where D <: Dict
  new_section("Material Models")
  models = Dict{String, ConstitutiveModels.ConstitutiveModel}()
  props  = Dict{String, Vector{Float64}}()
  states = Dict{String, Any}()
  for mat_name in keys(input_settings)
    @assert "model" in keys(input_settings[mat_name])
    @assert "properties" in keys(input_settings[mat_name])

    model_name = input_settings[mat_name]["model"]
    props_in   = input_settings[mat_name]["properties"]

    @info "Reading material $mat_name"
    @info "  Model      = $model_name"
    @info "  Properties = "
    for (key, val) in props_in
      @info "    $(rpad(key, 20)) = $val"
    end
    @info ""

    model, prop, state = eval(Meta.parse(model_name))(props_in)

    models[mat_name] = model
    props[mat_name]  = prop
    states[mat_name] = state
  end
  return models, props, states
end

function read_sections(input_settings, mesh, dof, models, props, states)
  new_section("Sections")
  block_ids = element_block_ids(mesh)
  sections = Dict()
  n = 1
  for section in input_settings
    @info "Reading Section $n"
    @warn "Defaulting to fully integrated element, e.g. q_degree = 2"
    @assert "block id" in keys(section)
    @assert "formulation" in keys(section)
    @assert "material" in keys(section)

    block_id    = section["block id"]
    formulation = section["formulation"]
    mat_name    = section["material"]
    q_degree    = 2 # TODO make this input somehow

    @assert block_id in block_ids
    @assert mat_name in keys(models)
    @assert mat_name in keys(props)
    @assert mat_name in keys(states)

    @info "  Block id      = $block_id"
    @info "  Formulation   = $formulation"
    @info "  Material name = $mat_name"
    @info ""

    conns     = element_connectivity(mesh, block_id)
    conns     = convert.(Int64, conns)
    conns     = Connectivity{size(conns, 1), size(conns, 2), Vector}(conns)
    
    elem_type = FiniteElementContainers.element_type(mesh, block_id)
    if formulation == "plane strain"
      form = FiniteElementContainers.PlaneStrain()
    elseif formulation == "three dimensional" 
      form = FiniteElementContainers.ThreeDimensional()
    else
      @assert false "Unsupported type"
    end

    # fspace  = NonAllocatedFunctionSpace(conns, dof, block_id, q_degree, elem_type)
    fspace  = NonAllocatedFunctionSpace(dof, conns, q_degree, elem_type)
    section = Section(fspace, form, models[mat_name], props[mat_name], states[mat_name]) 

    sections[Symbol("block_$block_id")] = section
    n = n + 1
  end

  # end_section("Sections")

  return NamedTuple(sections)
  # return NamedTuple{Tuple(keys(sections))}(values(sections))
end

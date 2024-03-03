abstract type Section{ID, FS, Form, Mod} end

function Base.size(section::Section)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NP = ConstitutiveModels.num_properties(section.model)
  NS = ConstitutiveModels.num_state_vars(section.model)
  return ND, NN, NP, NS
end

FiniteElementContainers.num_dimensions(sec::Section)           = FiniteElementContainers.num_dimensions(sec.fspace)
FiniteElementContainers.num_elements(sec::Section)             = FiniteElementContainers.num_elements(sec.fspace)
FiniteElementContainers.num_nodes_per_element(sec::Section)    = FiniteElementContainers.num_nodes_per_element(sec.fspace)
FiniteElementContainers.num_q_points(sec::Section)             = FiniteElementContainers.num_q_points(sec.fspace)

FiniteElementContainers.connectivity(sec::Section)             = connectivity(sec.fspace)
FiniteElementContainers.connectivity(sec::Section, e::Int)     = connectivity(sec.fspace, e)
FiniteElementContainers.dof_connectivity(sec::Section)         = dof_connectivity(sec.fspace)
FiniteElementContainers.dof_connectivity(sec::Section, e::Int) = dof_connectivity(sec.fspace, e)

quadrature_weights(section::Section) = ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe)
shape_function_gradients(section::Section) = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe)
shape_function_values(section::Section) = ReferenceFiniteElements.shape_function_values(section.fspace.ref_fe)

quadrature_weights(section::Section, q::Int) = ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe, q)
shape_function_gradients(section::Section, q::Int) = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
shape_function_values(section::Section, q::Int) = ReferenceFiniteElements.shape_function_values(section.fspace.ref_fe, q)


function interpolants(ref_fe::ReferenceFE, formulation, X_el, q::Int)
  ∇N_ξ  = ReferenceFiniteElements.shape_function_gradients(ref_fe, q)
  J     = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X  = (J_inv * ∇N_ξ')'
  JxW   = det(J) * ReferenceFiniteElements.quadrature_weights(ref_fe, q)
  G     = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return ∇N_X, JxW, G
end

# simple setup, eventually dispatch on setup type
function setup(section::Section, inputs::Dict)
  props_init, state_init = ConstitutiveModels.setup(section.model, inputs; type=SVector)

  NQ, NE = FiniteElementContainers.num_q_points(section), num_elements(section)

  # if length(state_init) == 0
  #   @warn "Assigning a single state variable to a path-indepent model"
  #   state_init = SVector{1, eltype(state_init)}(zero(eltype(state_init)))
  # end
  
  # scratch array for energy kernels
  Πs = ElementField{NQ, NE, Vector, Float64}(undef)
  Πs .= 0.0
  props = ElementField{length(props_init), NE, Matrix, eltype(props_init)}(undef)
  # state = QuadratureField{length(state_init), NQ, NE, StructVector, typeof(state_init)}(undef)

  # TODO add this type to FiniteElementContainers in some way as a wrapper
  state_old = Array{eltype(state_init), 3}(undef, length(state_init), NQ, NE)
  state_new = Array{eltype(state_init), 3}(undef, length(state_init), NQ, NE)

  for e in 1:NE
    props[:, e] = props_init
    for q in 1:NQ
      state_old[:, q, e] = state_init
      state_new[:, q, e] = state_init
    end
  end

  return props, state_old, state_new, Πs
end

include("TotalLagrangeSection.jl")

#######################################
# Some setup helpers

function setup_material(input_settings::D) where D <: Dict{Symbol, Any}
  model_name = input_settings[:model]
  model = eval(Meta.parse(model_name))()
  return model
end

function setup_sections(input_settings::D, mesh::FileMesh, dof) where D <: Vector{Dict{Symbol, Any}}
  new_section("Sections")
  block_ids = element_block_ids(mesh)
  sections = Dict{Symbol, Any}()
  state_old = Dict{Symbol, Any}()
  state_new = Dict{Symbol, Any}()
  props = Dict{Symbol, Any}()
  Πs = Dict{Symbol, Any}()

  n = 1

  for section in input_settings
    @info "Reading Section $n"
    @info "  type = $(section[:type])"
    @info "  formulation = $(section[:formulation])"
    @info "  material ="
    @info "    model name = $(section[:material][:model])"
    @info "    properties = "
    for (key, val) in section[:material][:properties]
      @info "      $(rpad(key, 20)) = $val"
    end
    @warn "Defaulting to fully integrated element, e.g. q_degree = 2"
    # q_degree = 1

    block_id = section[Symbol("block id")]
    @assert block_id in block_ids
    formulation = section[:formulation]
    mat_settings = section[:material]
    q_degree = 2 # TODO make this input somehow

    # setup up conn for fspace
    conns = convert.(Int64, element_connectivity(mesh, block_id))
    conns = Connectivity{size(conns), Vector}(conns)
    # conns = Connectivity{size(conns, 1), size(conns, 2), Vector, SVector}(conns)
    elem_type = FiniteElementContainers.element_type(mesh, block_id)
    fspace  = NonAllocatedFunctionSpace(dof, conns, q_degree, elem_type)

    # setup formulation for seeding/extracting material point stuff
    if formulation == "three dimensional"
      if FiniteElementContainers.num_dimensions(mesh) == 3
        form = FiniteElementContainers.ThreeDimensional()
      else
        @assert false "Unsupported formulation"
      end
    elseif formulation == "plane strain"
      form = FiniteElementContainers.PlaneStrain()
    else
      @assert false "Unsupported formulation type"
    end

    # setup material
    model = setup_material(mat_settings)
    
    # finally setup section
    section_type = eval(Meta.parse(section[:type]))
    section_temp = section_type(block_id, fspace, form, model)

    # TODO more to do here
    props_temp, state_old_temp, state_new_temp, Πs_temp = setup(section_temp, section[:material][:properties])

    # store section
    name = Symbol("section_$n")
    sections[name] = section_temp
    props[name] = props_temp
    state_old[name] = state_old_temp
    state_new[name] = state_new_temp
    Πs[name] = Πs_temp

    n = n + 1
  end
  return NamedTuple(sections), 
         ComponentArray(props), 
         ComponentArray(state_old), 
         ComponentArray(state_new), 
         ComponentArray(Πs)
end

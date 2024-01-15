abstract type Section{FS, Form, Mod} end

FiniteElementContainers.num_dimensions(sec::Section)           = FiniteElementContainers.num_dimensions(sec.fspace)
FiniteElementContainers.num_elements(sec::Section)             = FiniteElementContainers.num_elements(sec.fspace)
FiniteElementContainers.num_nodes_per_element(sec::Section)    = FiniteElementContainers.num_nodes_per_element(sec.fspace)
FiniteElementContainers.num_q_points(sec::Section)             = FiniteElementContainers.num_q_points(sec.fspace)

FiniteElementContainers.connectivity(sec::Section)             = connectivity(sec.fspace)
FiniteElementContainers.connectivity(sec::Section, e::Int)     = connectivity(sec.fspace, e)
FiniteElementContainers.dof_connectivity(sec::Section)         = dof_connectivity(sec.fspace)
FiniteElementContainers.dof_connectivity(sec::Section, e::Int) = dof_connectivity(sec.fspace, e)

function interpolants(ref_fe::ReferenceFE, formulation, X_el, q::Int)
  ∇N_ξ  = ReferenceFiniteElements.shape_function_gradients(ref_fe, q)
  J     = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X  = (J_inv * ∇N_ξ')'
  JxW   = det(J) * ReferenceFiniteElements.quadrature_weights(ref_fe, q)
  G     = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return ∇N_X, JxW, G
end

include("TotalLagrangeSection.jl")

#######################################
# Some setup helpers

function setup_material(input_settings::D) where D <: Dict{Symbol, Any}
  model_name = input_settings[:model]
  properties = input_settings[:properties]
  model, props, state = eval(Meta.parse(model_name))(properties)
  return model, props, state
end

function setup_sections(input_settings::D, mesh::FileMesh, dof) where D <: Vector{Dict{Symbol, Any}}
  new_section("Sections")
  block_ids = element_block_ids(mesh)
  sections = Dict{Symbol, Any}()
  n = 1

  for section in input_settings
    @show section
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
    block_id = section[Symbol("block id")]
    @assert block_id in block_ids
    formulation = section[:formulation]
    mat_settings = section[:material]
    q_degree = 2 # TODO make this input somehow

    # setup up conn for fspace
    conns = convert.(Int64, element_connectivity(mesh, block_id))
    conns = Connectivity{size(conns), Vector}(conns)
    elem_type = FiniteElementContainers.element_type(mesh, block_id)
    fspace  = NonAllocatedFunctionSpace(dof, conns, q_degree, elem_type)

    # setup formulation for seeding/extracting material point stuff
    if formulation == "plane strain"
      form = FiniteElementContainers.PlaneStrain()
    elseif formulation == "three dimensional" 
      form = FiniteElementContainers.ThreeDimensional()
    else
      @assert false "Unsupported type"
    end

    # setup material
    model, props, state = setup_material(mat_settings)

    # finally setup section
    section_type = eval(Meta.parse(section[:type]))
    section = section_type(fspace, form, model, props, state)

    # store section
    sections[Symbol("section_$n")] = section
    n = n + 1
  end
  return NamedTuple(sections)
end

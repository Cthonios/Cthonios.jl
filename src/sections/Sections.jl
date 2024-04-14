"""
$(TYPEDEF)
"""
abstract type Section{ID, FS, Form, Mod} end

"""
$(TYPEDSIGNATURES)
"""
function Base.size(section::Section)
  ND = FiniteElementContainers.num_dimensions(section)
  NN = FiniteElementContainers.num_nodes_per_element(section)
  NP = ConstitutiveModels.num_properties(section.model)
  NS = ConstitutiveModels.num_state_vars(section.model)
  return ND, NN, NP, NS
end

"""
$(TYPEDSIGNATURES)
Returns the dimensionality of the section, e.g. 2 or 3.
"""
FiniteElementContainers.num_dimensions(sec::Section)           = FiniteElementContainers.num_dimensions(sec.fspace)
"""
$(TYPEDSIGNATURES)
Returns the number of elements in the section.
"""
FiniteElementContainers.num_elements(sec::Section)             = FiniteElementContainers.num_elements(sec.fspace)
"""
$(TYPEDSIGNATURES)
Returns the number of nodes per element in the section.
"""
FiniteElementContainers.num_nodes_per_element(sec::Section)    = FiniteElementContainers.num_nodes_per_element(sec.fspace)
"""
$(TYPEDSIGNATURES)
Returns the number of quadrature points in the reference
element for a the section
"""
FiniteElementContainers.num_q_points(sec::Section)             = FiniteElementContainers.num_q_points(sec.fspace)
"""
$(TYPEDSIGNATURES)
Total element connectivity of the section.
"""
FiniteElementContainers.connectivity(sec::Section)             = connectivity(sec.fspace)
"""
$(TYPEDSIGNATURES)
Element connectivity for element ```e```
"""
FiniteElementContainers.connectivity(sec::Section, e::Int)     = connectivity(sec.fspace, e)
"""
$(TYPEDSIGNATURES)
Degree of freedom connectivity for all elements in the section.
"""
FiniteElementContainers.dof_connectivity(sec::Section)         = dof_connectivity(sec.fspace)
"""
$(TYPEDSIGNATURES)
Returns the degree of freedom connectivity.

For example if there is an element with connectivty [1, 2, 3, 4]
that is two dimensional, this method will return a connectivity of
[1, 2, 3, 4, 5, 6, 7, 8]
"""
FiniteElementContainers.dof_connectivity(sec::Section, e::Int) = dof_connectivity(sec.fspace, e)
"""
$(TYPEDSIGNATURES)
"""
quadrature_weights(section::Section) = ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe)
"""
$(TYPEDSIGNATURES)
"""
shape_function_gradients(section::Section) = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe)
"""
$(TYPEDSIGNATURES)
"""
shape_function_values(section::Section) = ReferenceFiniteElements.shape_function_values(section.fspace.ref_fe)
"""
$(TYPEDSIGNATURES)
"""
quadrature_weights(section::Section, q::Int) = ReferenceFiniteElements.quadrature_weights(section.fspace.ref_fe, q)
"""
$(TYPEDSIGNATURES)
"""
shape_function_gradients(section::Section, q::Int) = ReferenceFiniteElements.shape_function_gradients(section.fspace.ref_fe, q)
"""
$(TYPEDSIGNATURES)
"""
shape_function_values(section::Section, q::Int) = ReferenceFiniteElements.shape_function_values(section.fspace.ref_fe, q)

"""
$(TYPEDSIGNATURES)
Calculates interpolants on the fly given a reference element,
a formulation (e.g. ```ThreeDimensional```, ```PlaneStrain```, etc.),
coordinates for the element of interest in the reference (or current)
configuration and the quadrature point to evaluate at.
"""
function interpolants(ref_fe::ReferenceFE, formulation, X_el, q::Int)
  ∇N_ξ  = ReferenceFiniteElements.shape_function_gradients(ref_fe, q)
  J     = X_el * ∇N_ξ
  J_inv = inv(J)
  ∇N_X  = (J_inv * ∇N_ξ')'
  JxW   = det(J) * ReferenceFiniteElements.quadrature_weights(ref_fe, q)
  G     = FiniteElementContainers.discrete_gradient(formulation, ∇N_X)
  return ∇N_X, JxW, G
end

"""
sets up necessary section cache arrays for a single section given
constant properties
"""
function setup_section_cache(section, props_inputs)
  props_init, state_init = ConstitutiveModels.setup(section.model, props_inputs)
  NQ, NE = FiniteElementContainers.num_q_points(section), num_elements(section)

  Πs = ElementField{NQ, NE, Vector, Float64}(undef)
  Πs .= 0.0
  props = ElementField{length(props_init), NE, Matrix, eltype(props_init)}(undef)

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

"""
takes an unnamed but ordered set of properties and sections

Currently assumes all blocks in the mesh are being used

TODO add omit block type stuff
"""
function setup_section_caches(sections, props_inputs)
  @assert length(sections) == length(props_inputs) "Needs to be the same length"
  state_old = Dict{Symbol, Any}()
  state_new = Dict{Symbol, Any}()
  props = Dict{Symbol, Any}()
  Πs = Dict{Symbol, Any}()

  for (section, props_input) in zip(sections, props_inputs)
    key = Symbol("section_$(section.block_id)")
    p, s_old, s_new, pi = setup_section_cache(section, props_input)
    state_old[key] = s_old
    state_new[key] = s_new
    props[key] = p
    Πs[key] = pi
  end

  return ComponentArray(props), 
         ComponentArray(state_old), 
         ComponentArray(state_new), 
         ComponentArray(Πs)
end

include("TotalLagrangeSection.jl")

#######################################
# Some setup helpers
"""
$(TYPEDSIGNATURES)
"""
function setup_material(input_settings::D) where D <: Dict{Symbol, Any}
  model_name = input_settings[:model]
  model = eval(Meta.parse(model_name))(input_settings[:properties])
  return model
end

"""
$(TYPEDSIGNATURES)
"""
function setup_sections(input_settings::D, mesh::FileMesh, dof) where D <: Vector{Dict{Symbol, Any}}
  new_section("Sections")
  block_ids = element_block_ids(mesh)
  sections = Dict{Symbol, Any}()

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

    # store section
    name = Symbol("section_$n")
    sections[name] = section_temp

    n = n + 1
  end

  props_inputs = map(section -> section[:material][:properties], input_settings)
  props, state_old, state_new, Πs = setup_section_caches(values(sections), props_inputs)

  return NamedTuple(sections), 
         props,
         state_old,
         state_new,
         Πs
end
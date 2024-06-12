abstract type AbstractSection end
abstract type AbstractSectionInput end

struct TotalLagrangeSection <: AbstractSectionInput
  block_name::String
  formulation
  model
  q_degree::Int64
  props
end

struct TotalLagrangeSectionInternal{Fspace, Form, Model} <: AbstractSection
  # block_id::Int64
  block_name::String
  fspace::Fspace
  formulation::Form
  model::Model
end

function TotalLagrangeSectionInternal(
  mesh, dof, model, formulation, block_name, q_degree
)
  conns = convert.(Int64, element_connectivity(mesh, block_name))
  conns = Connectivity{size(conns), Vector}(conns)

  # get element type
  elem_type = FiniteElementContainers.element_type(mesh, block_name)

  # setup function space
  fspace = NonAllocatedFunctionSpace(dof, conns, q_degree, elem_type)

  section = TotalLagrangeSectionInternal{typeof(fspace), typeof(formulation), typeof(model)}(
    block_name, fspace, formulation, model
  )
  return section
end

"""
$(TYPEDSIGNATURES)
"""
function Base.size(section::AbstractSection)
  ND = FiniteElementContainers.num_dimensions(section.fspace)
  NN = FiniteElementContainers.num_nodes_per_element(section.fspace)
  NP = ConstitutiveModels.num_properties(section.model)
  NS = ConstitutiveModels.num_state_vars(section.model)
  return ND, NN, NP, NS
end

function modify_field_gradients(section, ∇u)
  return FiniteElementContainers.modify_field_gradients(
    section.formulation, ∇u
  )
end

# Sections
# TODO make more amenable to multi dispatch
function setup_section_cache(section_in, section::TotalLagrangeSectionInternal)
  props_inputs = section_in.props
  props_init, state_init = ConstitutiveModels.setup(section.model, props_inputs)
  NQ, NE = FiniteElementContainers.num_q_points(section.fspace), num_elements(section.fspace)

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

function setup_section_caches(sections_in, sections)
  # @assert length(sections) == length(props_inputs) "Needs to be the same length"
  state_old = Dict{Symbol, Any}()
  state_new = Dict{Symbol, Any}()
  props = Dict{Symbol, Any}()
  Πs = Dict{Symbol, Any}()

  # for (section, props_input) in zip(sections, props_inputs)
  for (section_in, section) in zip(sections_in, sections)
    key = Symbol("section_$(section.block_name)")
    p, s_old, s_new, pi = setup_section_cache(section_in, section)
    state_old[key] = s_old
    state_new[key] = s_new
    props[key] = p
    Πs[key] = pi
  end

  return ComponentArray(Πs),
         ComponentArray(props),
         ComponentArray(state_old),
         ComponentArray(state_new)
end

function setup_sections(mesh, dof, sections)
  sections_out = Dict{Symbol, Any}()
  for section in sections
    section_name = Symbol("section_$(section.block_name)")
    sections_out[section_name] = TotalLagrangeSectionInternal(
      mesh, dof, 
      section.model, section.formulation, 
      section.block_name, section.q_degree
    )
  end
  return NamedTuple(sections_out)
end

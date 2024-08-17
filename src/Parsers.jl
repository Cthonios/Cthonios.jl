required_keys = Symbol[
  Symbol("dirichlet boundary conditions"),
  Symbol("neumann boundary conditions"),
  # :mesh,
  :sections,
  :solver,
  Symbol("time stepper"),
  #
  :Domain,
  :Objective,
  :PostProcessor
]

struct InputFileSectionNotFound <: Exception
  required_key::Symbol
end

function Base.show(io::IO, e::InputFileSectionNotFound)
  println(io, "$(e.required_key) not found in input file.")
end

function check_keys(inputs::Dict{Symbol, Any})
  for required_key in required_keys
    if required_key ∉ keys(inputs)
      throw(InputFileSectionNotFound(required_key))
    end
  end
  return nothing
end

function get_dict(inputs::Dict{Symbol, Any}, key::Symbol)::Dict{Symbol, Any}
  inputs[key]
end

function get_vector_dict(inputs::Dict{Symbol, Any}, key::Symbol)::Vector{Dict{Symbol, Any}}
  inputs[key]
end

function load_inputs(input_file::String)::Dict{Symbol, Any}
  YAML.load_file(input_file; dicttype=Dict{Symbol, Any})
end

function setup_group(item::Dict{Symbol, Any})
  types = keys(item) |> collect
  @assert length(types) == 1
  type = types[1]
  inputs = item[type]
  eval(type)(inputs)
end

function setup_group(group::Vector{T}, type) where T
  new_items = type[]
  for item in group
    types = keys(item) |> collect
    @assert length(types) == 1
    type = types[1]
    inputs = item[type]
    new_item = eval(type)(inputs)
    push!(new_items, new_item)
  end
  new_items
end

function parse_input_file(input_file::String)
  inputs = load_inputs(input_file)
  check_keys(inputs)
  # read inputs
  dbc_inputs = get_vector_dict(inputs, Symbol("dirichlet boundary conditions"))
  nbc_inputs = get_vector_dict(inputs, Symbol("neumann boundary conditions"))

  # mesh_inputs = get_dict(inputs, :mesh)
  section_inputs = get_vector_dict(inputs, :sections)
  solver_inputs = get_dict(inputs, :solver)
  time_stepper_inputs = get_dict(inputs, Symbol("time stepper"))
  # setup groups based on inputs
  dbcs = setup_group(dbc_inputs, DirichletBC)
  nbcs = setup_group(nbc_inputs, NeumnanBC)
  # mesh = setup_group(mesh_inputs)
  sections = setup_group(section_inputs, TotalLagrangeSection)
  solver = setup_group(solver_inputs)
  time_stepper = setup_group(time_stepper_inputs)
  #
  # domain = setup_group(domain_inputs)
  domain = Domain(
    inputs[:Domain], time_stepper,
    sections, dbcs, nbcs,
    solver, 2 # TODO 2 is for 2 dofs need to change this
  )
  pp = PostProcessor(inputs[:PostProcessor], domain)
  obj = Objective(inputs[:Objective], domain, pp)
end
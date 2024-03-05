abstract type AbstractDomainCache end

abstract type AbstractDomain{
  Dof,
  Funcs,
  BCNodes,
  BCDofs,
  BCFuncIDs,
  Sections
} end

include("QuasiStaticDomain.jl")

# Parsing
get_domain_bcs_inputs(input_settings)::Dict{Symbol, Any} = input_settings[:boundary_conditions]
get_domain_displacement_bcs_inputs(input_settings)::Vector{Dict{Symbol, Any}} = get_domain_bcs_inputs(input_settings)[:displacement] 
get_domain_sections_inputs(input_settings)::Vector{Dict{Symbol, Any}} = input_settings[:sections]
get_domain_time_stepper_inputs(input_settings)::Dict{Symbol, Any} = input_settings[:time_stepper]

required_domain_keys = String[
  "mesh",
  "functions",
  "boundary conditions",
  "materials",
  "sections"
]

create_mesh(type, input_settings) = FileMesh(type, input_settings)

function read_mesh(input_settings)#::FileMesh{ExodusDatabase{Int32, Int32, Int32, Float64}} where D <: Dict
  @assert isfile(input_settings[:mesh][Symbol("file name")])
  @info "Mesh file = $(input_settings[:mesh][Symbol("file name")])"
  @info "Opening mesh"
  type = eval(Meta.parse(input_settings[:mesh][:type]))
  mesh = create_mesh(type, input_settings[:mesh][Symbol("file name")])
  return mesh
end

function read_coordinates(mesh)#::Matrix{Float64}
  @info "Reading coordinates"
  coords = FiniteElementContainers.coordinates(mesh)
  coords = NodalField{size(coords), Vector}(coords)
  @info "Read coordinates into field of type $(typeof(coords))"
  return coords
end

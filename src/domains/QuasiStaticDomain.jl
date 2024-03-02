struct QuasiStaticDomainCache{V1, V2, V3, V4, V5, V6, V7, V8} <: AbstractDomainCache
  X::V1
  U::V2
  state::V3
  props::V4
  Π::V5
  Πs::V6 # for energy kernels
  f::V7  # for internal force
  V::V8 # for working with krylov methods
end

function Base.similar(cache::QuasiStaticDomainCache)
  @unpack X, U, state, props, Π, Πs, f, V = cache
  return QuasiStaticDomainCache(
    similar(X), similar(U), similar(state),
    similar(props), similar(Π), similar(Πs), similar(f), similar(V)
  )
end

function unpack(cache::QuasiStaticDomainCache)
  @unpack X, U, state, props, Π, Πs, f, V = cache
  return X, U, state, props, Π, Πs, f, V
end

struct QuasiStaticDomain{
  Dof,
  Funcs,
  BCNodes,
  BCDofs,
  BCFuncIDs,
  Sections,
  Time,
  DomainCache
} <: AbstractDomain{Dof, Funcs, BCNodes, BCDofs, BCFuncIDs, Sections, Time, DomainCache}
  dof::Dof
  funcs::Funcs
  bc_nodes::BCNodes
  bc_dofs::BCDofs
  bc_func_ids::BCFuncIDs
  sections::Sections
  time::Time
  domain_cache::DomainCache
end

function QuasiStaticDomain(input_settings::D) where D <: Dict{Symbol, Any}

  # global stuff
  # gather all functions from bcs
  funcs = setup_functions(map(x -> x[:function], get_domain_displacement_bcs_inputs(input_settings)))

  # read mesh and set up some fields
  new_section("Mesh")
  mesh_file = read_mesh(input_settings)
  ND        = FiniteElementContainers.num_dimensions(mesh_file) |> Int64
  NNodes    = FiniteElementContainers.num_nodes(mesh_file) |> Int64
  coords    = read_coordinates(mesh_file)
  dof       = DofManager{ND, NNodes, Vector{Float64}}()
  
  # bc setup
  disp_bc_nodes, disp_bc_dofs, disp_bc_func_ids = setup_displacement_bcs(
    get_domain_displacement_bcs_inputs(input_settings), mesh_file, ND
  )
  
  # section setup
  sections, props, state, Πs = setup_sections(get_domain_sections_inputs(input_settings), mesh_file, dof)

  # time stepper setup
  time = ConstantTimeStepper(get_domain_time_stepper_inputs(input_settings))

  # cache setup
  # X = copy(coords)
  U = FiniteElementContainers.create_fields(dof)
  f = FiniteElementContainers.create_fields(dof).vals
  V = FiniteElementContainers.create_fields(dof)
  domain_cache = QuasiStaticDomainCache(coords, U, state, props, zeros(Float64, 1), Πs, f, V)

  return QuasiStaticDomain(
    dof, funcs, 
    disp_bc_nodes, disp_bc_dofs, disp_bc_func_ids,
    sections, time, domain_cache
  )
end

QuasiStaticDomain(input_file::String, key) =
QuasiStaticDomain(parse_input_file(input_file)[:domains][key])

function Base.show(io::IO, domain::QuasiStaticDomain)
  print(io, "\n    QuasiStaticDomain\n")
  print(io, "      Sections\n")
  for (key, val) in pairs(domain.sections)
    println(io, "        $key")
    println(io, val)
  end
  print(io, "      TimeStepper\n")
  print(io, "        Type = $(typeof(domain.time).name.name)")
end

# access methods
coordinates(d::QuasiStaticDomain)  = d.coords
dof_manager(d::QuasiStaticDomain)  = d.dof
sections(d::QuasiStaticDomain)     = d.sections
time_stepper(d::QuasiStaticDomain) = d.time

# FEM container methods
# DofManager hooks
create_fields(d::QuasiStaticDomain) = FiniteElementContainers.create_fields(d.dof)
create_unknowns(d::QuasiStaticDomain) = FiniteElementContainers.create_unknowns(d.dof)

"""
This methods assumes the bc dofs and func ids are already
properly set in the domain coming in
"""
function update_unknown_dofs!(d::QuasiStaticDomain)
  # update the dofs
  FiniteElementContainers.update_unknown_dofs!(d.dof, d.bc_dofs)
end

function update_bcs!(U, domain::QuasiStaticDomain, Xs)
  for (node, dof, func_id) in zip(domain.bc_nodes, domain.bc_dofs, domain.bc_func_ids)
    X = @views Xs[:, node]
    t = domain.time.current_time
    U[dof] = domain.funcs[func_id](X, t)
  end
end

function update_unknowns!(U, domain::QuasiStaticDomain, Uu)
  FiniteElementContainers.update_fields!(U, domain.dof, Uu)
end

function update_fields!(U, domain::QuasiStaticDomain, Xs, Uu)
  update_bcs!(U, domain, Xs)
  update_unknowns!(U, domain, Uu)
end

# Time stepping hooks
step!(domain::QuasiStaticDomain) = step!(domain.time)

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

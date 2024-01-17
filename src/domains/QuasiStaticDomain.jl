# struct QuasiStaticDomain{
#   Coords    <: NodalField,
#   Dof       <: DofManager,
#   # Funcs     <: ScalarFunctionWrapperArray,
#   Funcs, # type eventually, above is currently give instability
#   BCNodes   <: AbstractArray{<:Integer, 1},
#   BCDofs    <: AbstractArray{<:Integer, 1},
#   BCFuncIDs <: AbstractArray{<:Integer, 1},
#   Sections  <: NamedTuple,
#   Assembler <: StaticAssembler,
#   Time      <: TimeStepper
# } <: Domain
struct QuasiStaticDomain{
  Coords,
  Dof,
  Funcs,
  BCNodes,
  BCDofs,
  BCFuncIDs,
  Sections,
  Assembler,
  Time
}
  coords::Coords
  dof::Dof
  funcs::Funcs
  bc_nodes::BCNodes
  bc_dofs::BCDofs
  bc_func_ids::BCFuncIDs
  sections::Sections
  assembler::Assembler
  time::Time
end

get_domain_bcs_inputs(input_settings)::Dict{Symbol, Any} = input_settings[:boundary_conditions]
get_domain_displacement_bcs_inputs(input_settings)::Vector{Dict{Symbol, Any}} = get_domain_bcs_inputs(input_settings)[:displacement] 
get_domain_sections_inputs(input_settings)::Vector{Dict{Symbol, Any}} = input_settings[:sections]
get_domain_time_stepper_inputs(input_settings)::Dict{Symbol, Any} = input_settings[:time_stepper]

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
  @info "Setting up vectorized DofManager"
  dof       = DofManager{ND, NNodes, Vector{Float64}}()
  
  # bc setup
  disp_bc_nodes, disp_bc_dofs, disp_bc_func_ids = setup_displacement_bcs(
    get_domain_displacement_bcs_inputs(input_settings), mesh_file, ND
  )
  
  # section setup
  sections = setup_sections(get_domain_sections_inputs(input_settings), mesh_file, dof)

  # assembler setup
  assembler = StaticAssembler(dof, map(x -> x.fspace, values(sections)))

  # time stepper setup
  time = ConstantTimeStepper(get_domain_time_stepper_inputs(input_settings))

  return QuasiStaticDomain(
    coords, dof, funcs, 
    disp_bc_nodes, disp_bc_dofs, disp_bc_func_ids,
    sections, assembler, time
  )
end

function QuasiStaticDomain_old(input_settings::D) where D <: Dict{Symbol, Any}

  # read mesh and do some setup
  mesh_file = read_mesh(input_settings["mesh"])
  ND        = FiniteElementContainers.num_dimensions(mesh_file)
  coords    = read_coordinates(mesh_file)
  NFields   = FiniteElementContainers.num_dimensions(mesh_file) |> Int64
  NNodes    = size(coords, 2)
  dof       = DofManager{NFields, NNodes, Vector}() # TODO issue with Vector here
  funcs     = read_functions(input_settings["functions"], ND)
  # add defaults to funcs
  # TODO add unitful stuff eventually
  # TODO figure out how to add defaults
  # zero_func = ScalarFunction{2, Float64, Float64, Float64}((x, t) -> 0.0)
  # funcs = (; funcs..., :zero => zero_func)
  # read boundary conditions with mesh and funcs as inputs
  bcs = read_bcs(input_settings["boundary conditions"], mesh_file, funcs, NFields)
  # read material models and set them up, then set up sections, and finally the assembler
  models, props, states = read_materials(input_settings["materials"])
  sections              = read_sections(input_settings["sections"], mesh_file, dof, models, props, states)  
  assembler             = StaticAssembler(dof, map(x -> x.fspace, values(sections)))

  # TODO refine this to take in other time steppers
  time = ConstantTimeStepper(input_settings["time stepper"])

  # postprocess funcs into a function wrapper array
  # funcs = values(funcs) |> collect |> ScalarFunctionWrapperArray
  # funcs = ntuple(funcs...)
  # get all bc dofs and bc func ids
  # TODO need to do check on conflicting boundary conditions
  bc_nodes = mapreduce(x -> x.nodes, vcat, bcs)
  bc_dofs = mapreduce(x -> x.dofs, vcat, bcs)
  bc_func_ids = mapreduce(x -> fill(x.func_id, length(x.dofs)), vcat, bcs)

  # filter out repeat dofs, this is a first come first serve
  unique_bc_dof_ids = unique(i -> bc_dofs[i], eachindex(bc_dofs))

  unique!(bc_dofs)
  # TODO maybe add sorting to the node ids for faster access?
  bc_nodes = bc_nodes[unique_bc_dof_ids]
  bc_func_ids = bc_func_ids[unique_bc_dof_ids]

  # now sort them
  # sorted_indices = sortperm(bc_dofs)
  # bc_dofs = 
  return QuasiStaticDomain(coords, dof, funcs, bc_nodes, bc_dofs, bc_func_ids, sections, assembler, time)
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
FiniteElementContainers.create_fields(d::QuasiStaticDomain) = create_fields(d.dof)
FiniteElementContainers.create_unknowns(d::QuasiStaticDomain) = create_unknowns(d.dof)

"""
This methods assumes the bc dofs and func ids are already
properly set in the domain coming in
"""
function FiniteElementContainers.update_unknown_dofs!(d::QuasiStaticDomain)
  update_unknown_dofs!(d.dof, d.bc_dofs)
  update_unknown_dofs!(d.assembler, map(x -> x.fspace, d.sections), d.bc_dofs)
end

step!(domain::QuasiStaticDomain) = step!(domain.time)

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

function FiniteElementContainers.update_fields!(U, domain::QuasiStaticDomain, Xs, Uu)
  update_bcs!(U, domain, Xs)
  update_unknowns!(U, domain, Uu)
end

"""
Energy method that is a wrapper for nodal fields
"""
function energy(domain::QuasiStaticDomain, U::V1, X::V2) where {V1 <: NodalField, V2}
  W = 0.0 # Unitful error here # TODO 
  for section in domain.sections
    W = W + energy(section, U, X)
  end
  W
end

"""
Method to be called by a user
"""
function energy(domain::QuasiStaticDomain, Uu::V1, X::V2) where {V1 <: AbstractVector, V2}
  U = create_fields(domain)
  update_fields!(U, domain, X, Uu)
  energy(domain, U, X)
end

"""
In place residual method
"""
function residual!(domain::QuasiStaticDomain, Uu::V1, X::V2, U) where {V1 <: AbstractArray, V2}
  domain.assembler.residuals .= 0.0
  update_unknowns!(U, domain, Uu)
  for section in domain.sections
    residual!(domain.assembler, section, U, X)
  end
end

"""
In place hvp method

TODO right now this abuses the memory for residual. Maybe this is ok
"""
function hvp!(domain::QuasiStaticDomain, Uu::V1, Vv::V2, X::V3, U, V) where {V1 <: AbstractArray, V2, V3}
  domain.assembler.residuals .= 0.0
  update_unknowns!(U, domain, Uu)
  # update_unknowns!(V, domain, Vv)
  update_unknowns!(V, domain, Vv)
  for section in domain.sections
    hvp!(domain.assembler, section, U, V, X)
  end
end

"""
In place stiffness method
"""
function stiffness!(domain::QuasiStaticDomain, Uu::V1, X::V2, U) where {V1 <: AbstractArray, V2}
  domain.assembler.stiffnesses .= 0.0
  update_unknowns!(U, domain, Uu)
  for (n, section) in enumerate(domain.sections)
    stiffness!(domain.assembler, section, U, X, n)
  end
end

"""
Out of place residual method.
"""
function residual(domain::QuasiStaticDomain, Uu::V1, X::V2) where {V1 <: AbstractArray, V2}
  U = create_fields(domain)
  update_fields!(U, domain, X, Uu)
  residual!(domain, Uu, X, U)
  return domain.assembler.residuals[domain.dof.unknown_dofs]
end

"""
Out of place hvp method.
"""
function hvp(domain::QuasiStaticDomain, Uu::V1, Vv::V2, X::V3) where {V1 <: AbstractArray, V2, V3}
  U = create_fields(domain)
  V = create_fields(domain)
  update_fields!(U, domain, X, Uu)
  # update_unknowns!(V, domain, Vv)
  update_fields!(V, domain, X, Vv)
  hvp!(domain, Uu, Vv, X, U, V)
  return domain.assembler.residuals[domain.dof.unknown_dofs]
end

"""
Out of place stiffness method.
"""
function stiffness(domain::QuasiStaticDomain, Uu::V1, X::V2) where {V1 <: AbstractArray, V2}
  U = create_fields(domain)
  update_fields!(U, domain, X, Uu)
  stiffness!(domain, Uu, X, U)
  # return sparse(domain.assembler) |> symmetric
  K = sparse(domain.assembler)
  return 0.5 * (K + K')
end

# helper methods when we don't care about shape optimization
energy(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector =
energy(domain, Uu, domain.coords)
residual(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector =
residual(domain, Uu, domain.coords)
hvp(domain::QuasiStaticDomain, Uu::V1, Vv::V2) where {V1 <: AbstractVector, V2 <: AbstractVector} =
hvp(domain, Uu, Vv, domain.coords)
stiffness(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector =
stiffness(domain, Uu, domain.coords)
##########################################################################
# Parsing and setup helpers

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

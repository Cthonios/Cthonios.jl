abstract type Domain <: AbstractCthoniosType end

struct QuasiStaticDomain{
  Coords    <: NodalField,
  Dof       <: DofManager,
  Sections  <: NamedTuple,
  Assembler <: StaticAssembler,
  BCs       <: NamedTuple,
  Time      <: TimeStepper,
  Post      <: PostProcessor
} <: Domain
  coords::Coords
  dof::Dof
  sections::Sections
  assembler::Assembler
  bcs::BCs
  time::Time
  post_processor::Post
end

function QuasiStaticDomain(input_settings::D) where D <: Dict
  # TODO make variable maybe?
  # arr_type = Vector
    
  for key in required_domain_keys
    @assert key in keys(input_settings)
  end

  # break up input settings
  mesh_settings   = input_settings["mesh"]
  bc_settings     = input_settings["boundary conditions"]
  mat_settings    = input_settings["materials"]
  sec_settings    = input_settings["sections"]
  time_settings   = input_settings["time stepper"]

  # set stuff up
  mesh_file             = read_mesh(mesh_settings)
  coords                = read_coordinates(mesh_file)
  NFields               = FiniteElementContainers.num_dimensions(mesh_file) |> Int64
  NNodes                = size(coords, 2)
  dof                   = DofManager{NFields, NNodes, Vector}() # TODO issue with Vector here
  bcs                   = read_bcs(bc_settings, mesh_file)
  models, props, states = read_materials(mat_settings)
  sections              = read_sections(sec_settings, mesh_file, dof, models, props, states)  
  assembler             = StaticAssembler(dof, map(x -> x.fspace, values(sections)))
  time                  = ConstantTimeStepper(time_settings)
  post_processor        = PostProcessor(mesh_file, "output.e", dof, FiniteElementContainers.num_dimensions(mesh_file) |> Int64) # TODO change this

  return QuasiStaticDomain(coords, dof, sections, assembler, bcs, time, post_processor)
end

function Base.show(io::IO, domain::QuasiStaticDomain)
  print(io, "\n    QuasiStaticDomain\n")
  print(io, "      Sections\n")
  for (key, val) in pairs(domain.sections)
    println(io, "        $key")
    println(io, val)
  end
  println(io, "      Boundary conditions")
  for (key, val) in pairs(domain.bcs)
    println(io, "        $key")
    println(io, val)
  end
  println(io, "      PostProcessor")
  println(io, domain.post_processor)
end

bounbdary_conditions(d::QuasiStaticDomain) = d.bcs
coordinates(d::QuasiStaticDomain)          = d.coords
dof_manager(d::QuasiStaticDomain)          = d.dof
sections(d::QuasiStaticDomain)             = d.sections

num_dimensions(d::QuasiStaticDomain) = FiniteElementContainers.num_fields(coordinates(d))

create_fields(d::QuasiStaticDomain)         = FiniteElementContainers.create_fields(d.dof)
create_fields(d::QuasiStaticDomain, type)   = FiniteElementContainers.create_fields(d.dof, type)
create_unknowns(d::QuasiStaticDomain)       = FiniteElementContainers.create_unknowns(d.dof)
create_unknowns(d::QuasiStaticDomain, type) = FiniteElementContainers.create_unknowns(d.dof, type)

# TODO place where units won't work
function update_unknown_ids!(dof::DofManager, bcs)
  for bc in bcs
    FiniteElementContainers.update_unknown_ids!(dof, bc.nodes, bc.dof)
  end
  resize!(dof.unknown_indices, sum(dof.is_unknown))
  dof.unknown_indices .= FiniteElementContainers.dof_ids(dof)[dof.is_unknown]
end

function update_unknown_ids!(d::QuasiStaticDomain)

  dof = dof_manager(d)
  for bc in bounbdary_conditions(d)
    FiniteElementContainers.update_unknown_ids!(dof_manager(d), bc.nodes, bc.dof)
  end
  
  resize!(dof.unknown_indices, sum(dof.is_unknown))
  dof.unknown_indices .= FiniteElementContainers.dof_ids(dof)[dof.is_unknown]
end

function update_bcs!(U, coords, time, bcs)
  for bc in bcs
    for node in bc.nodes
      # get coords
      X = @views coords[:, node]
      t = time.current_time
      U[bc.dof, node] = bc.func(X, t)
    end
  end
end

function update_bcs!(U, d::QuasiStaticDomain)
  for bc in d.bcs
    for node in bc.nodes
      # get coords
      X = @views d.coords[:, node]
      t = d.time.current_time
      U[bc.dof, node] = bc.func(X, t)
    end
  end
end

# for structural optimization problems
function update_bcs!(U, X, d::QuasiStaticDomain)
  for bc in d.bcs
    for node in bc.nodes
      # get coords
      x = @views X[:, node]
      t = d.time.current_time
      U[bc.dof, node] = bc.func(x, t)
    end
  end
end

# update_fields!(domain::QuasiStaticDomain) = FiniteElementContainers.update_fields!(domain.U, domain.dof, domain.Uu)
update_fields!(U::V1, domain::QuasiStaticDomain, Uu::V2) where {V1 <: NodalField, V2 <: AbstractArray{<:Number, 1}} =
FiniteElementContainers.update_fields!(U, domain.dof, Uu)
update_unknowns!(U::V1, domain::QuasiStaticDomain, Uu::V2) where {V1 <: NodalField, V2 <: AbstractArray{<:Number, 1}} =
FiniteElementContainers.update_fields!(U, domain.dof, Uu)
################################################################

# useful for adjoint calculations and other AD
function energy(domain::QuasiStaticDomain, X::V1, U::V2) where {V1 <: Union{Matrix, NodalField}, V2 <: NodalField}
  W = 0.0 # Unitful error here # TODO 
  for section in domain.sections
    W = W + energy(section, X, U)
  end
  W
end

function energy(domain, X::V1, Uu::V2) where {V1 <: Union{Matrix, NodalField}, V2 <: AbstractVector}
  U = create_fields(domain, eltype(Uu))
  update_bcs!(U, domain)
  update_unknowns!(U, domain, Uu)
  energy(domain, X, U)
end

energy(domain, Uu) = energy(domain, domain.coords, Uu)

function residual(domain::QuasiStaticDomain, X::V1, U::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}
  R = Cthonios.create_fields(domain, eltype(U))
  for section in domain.sections
    residual!(R, section, X, U)
  end
  return R[domain.dof.is_unknown]
end

function residual(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
  U = Cthonios.create_fields(domain, eltype(Uu))
  update_fields!(U, domain, Uu)
  update_bcs!(U, domain.coords, domain.time, domain.bcs)
  return residual(domain, domain.coords, U)
end

function stiffness(domain::QuasiStaticDomain, X::V, U::V) where V <: NodalField
  domain.assembler.K .= 0.0 # Unitful issue here
  for section in domain.sections
    stiffness!(domain.assembler.K, section, X, U)
  end
  K = domain.assembler.K[domain.dof.is_unknown, domain.dof.is_unknown]
  return 0.5 * (K + K')
  # return K
end

function stiffness(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector{<:Number}
  U = Cthonios.create_fields(domain, eltype(Uu))
  update_fields!(U, domain, Uu)
  update_bcs!(U, domain.coords, domain.time, domain.bcs)
  return stiffness(domain, domain.coords, U)
end

# AD methods
grad_energy_x(backend, domain::QuasiStaticDomain, x, u) = AD.gradient(backend, z -> energy(domain, z, u), x)[1]
grad_energy_u(backend, domain::QuasiStaticDomain, x, u) = AD.gradient(backend, z -> energy(domain, x, z), u)[1]

# note these return a method that takes in a direction v
hvp_energy_x(backend, domain::QuasiStaticDomain, x, u) = 
AD.pushforward_function(backend, z -> grad_energy_x(backend, domain, z, u), x)
hvp_energy_u(backend, domain::QuasiStaticDomain, x, u) = 
AD.pushforward_function(backend, z -> grad_energy_u(backend, domain, x, z), u)

##########################################################################
# Parsing and setup helpers

required_domain_keys = String[
  "mesh",
  "boundary conditions",
  "materials",
  "sections"
]

function read_mesh(input_settings::D) where D <: Dict
  @assert "file name" in keys(input_settings)
  @assert isfile(input_settings["file name"])
  @info "Mesh file = $(input_settings["file name"])"
  @info "Opening mesh"
  mesh = FileMesh(ExodusDatabase, input_settings["file name"])
  @info @show mesh.mesh_obj
  return mesh
end

function read_coordinates(mesh)
  @info "Reading coordinates"
  coords = FiniteElementContainers.coordinates(mesh)
  coords = NodalField{size(coords, 1), size(coords, 2), Vector}(coords)
  @info "Read coordinates into field of type $(typeof(coords))"
  return coords
end


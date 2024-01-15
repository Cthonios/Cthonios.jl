abstract type Domain <: AbstractCthoniosType end

function energy end
function energy_gradient_u end
function energy_gradient_x end
function energy_hvp_u end

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


create_fields(d::QuasiStaticDomain, ::KernelAbstractionsBackend{Nothing}) = 
FiniteElementContainers.create_fields(d.dof)
create_fields(d::QuasiStaticDomain, type, ::KernelAbstractionsBackend{Nothing}) = 
FiniteElementContainers.create_fields(d.dof, type)
create_unknowns(d::QuasiStaticDomain, ::KernelAbstractionsBackend{Nothing}) = 
FiniteElementContainers.create_unknowns(d.dof)
create_unknowns(d::QuasiStaticDomain, type, ::KernelAbstractionsBackend{Nothing}) = 
FiniteElementContainers.create_unknowns(d.dof, type)


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
  
  # TODO add update_unknown_ids! for assembler here

  resize!(dof.unknown_indices, sum(dof.is_unknown))
  dof.unknown_indices .= FiniteElementContainers.dof_ids(dof)[dof.is_unknown]
end

function update_bcs!(U, coords, t, bcs, ::KernelAbstractionsBackend{Nothing})
  for bc in bcs
    for node in bc.nodes
      # get coords
      X = @views coords[:, node]
      U[bc.dof, node] = bc.func(X, t)
    end
  end
end

# update_fields!(domain::QuasiStaticDomain) = FiniteElementContainers.update_fields!(domain.U, domain.dof, domain.Uu)
update_fields!(U::V1, domain::QuasiStaticDomain, Uu::V2, ::KernelAbstractionsBackend{Nothing}) where {V1 <: NodalField, V2 <: AbstractArray{<:Number, 1}} =
FiniteElementContainers.update_fields!(U, domain.dof, Uu)
update_unknowns!(U::V1, domain::QuasiStaticDomain, Uu::V2, ::KernelAbstractionsBackend{Nothing}) where {V1 <: NodalField, V2 <: AbstractArray{<:Number, 1}} =
FiniteElementContainers.update_fields!(U, domain.dof, Uu)
################################################################

# useful for adjoint calculations and other AD
# function energy(domain::QuasiStaticDomain, X::V1, U::V2) where {V1 <: Union{Matrix, NodalField}, V2 <: NodalField}
# TODO wrap for shpa optimization adjoint later
# TODO also add state var stuff
function energy(domain::QuasiStaticDomain, U::V1, X::V2, ka_backend::KernelAbstractionsBackend) where {V1 <: NodalField, V2}
  W = 0.0 # Unitful error here # TODO 
  for section in domain.sections
    W = W + energy(section, U, X, ka_backend)
  end
  W
end

function energy(domain::QuasiStaticDomain, Uu::V, p, t, ka_backend::KernelAbstractionsBackend) where V <: AbstractVector
  U = create_fields(domain, eltype(Uu))
  update_bcs!(U, domain.coords, t, domain.bcs, ka_backend)
  update_unknowns!(U, domain, Uu, ka_backend)
  energy(domain, U, domain.coords, ka_backend)
end

function energy(domain::QuasiStaticDomain, Uu::V1, X::V2, ka_backend::KernelAbstractionsBackend) where {V1 <: AbstractVector, V2}
  U = create_fields(domain, eltype(Uu), ka_backend)
  update_bcs!(U, domain.coords, domain.time.current_time, domain.bcs, ka_backend)
  update_unknowns!(U, domain, Uu, ka_backend)
  energy(domain, U, X, ka_backend)
end

energy(domain::QuasiStaticDomain, Uu::V, ka_backend::KernelAbstractionsBackend) where V <: AbstractVector =
energy(domain, Uu, domain.coords, ka_backend)

function energy!(W::V1, domain::QuasiStaticDomain, Uu::V1, U::V2, ka_backend::KernelAbstractionsBackend) where {V1, V2}
  update_bcs!(U, domain.coords, domain.time.current_time, domain.bcs, ka_backend)
  update_unknowns!(U, domain, Uu, ka_backend)
  W[1] = energy(domain, U, domain.coords, ka_backend)
  nothing
end

function residual(
  domain::QuasiStaticDomain, U::V,
  ka_backend::KernelAbstractionsBackend
) where V <: AbstractArray
  # R = Cthonios.create_fields(domain, eltype(U))
  for section in domain.sections
    residual!(domain.assembler, section, U, domain.coords, ka_backend)
  end
  return R[domain.dof.is_unknown]
end

"""
This method assumes U comes in with BCs set
"""
function residual!(R::V1, domain::QuasiStaticDomain, Uu::V2, U::V1, ka_backend::KernelAbstractionsBackend) where {V1 <: NodalField, V2 <: AbstractVector}
  update_fields!(U, domain, Uu, ka_backend)
  for section in domain.sections
    residual!(R, section, U, domain.coords, ka_backend)
  end
end

function residual(
  domain::QuasiStaticDomain, Uu::V,
  ka_backend::KernelAbstractionsBackend
) where V <: AbstractVector
  U = Cthonios.create_fields(domain, eltype(Uu))
  update_fields!(U, domain, Uu, ka_backend)
  update_bcs!(U, domain.coords, domain.time.current_time, domain.bcs, ka_backend)
  return residual(domain, U, ka_backend)
end

function stiffness!(K, domain::QuasiStaticDomain, U::V, ka_backend::KernelAbstractionsBackend) where V <: NodalField
  K .= 0.0 # Unitful issue here
  for section in domain.sections
    stiffness!(K, section, U, domain.coords, ka_backend)
  end
end

function stiffness(domain::QuasiStaticDomain, U::V, ka_backend::KernelAbstractionsBackend) where V <: NodalField
  stiffness!(domain.assembler.K, domain, U, ka_backend)
  K = domain.assembler.K[domain.dof.is_unknown, domain.dof.is_unknown]
  return 0.5 * (K + K')
end

function stiffness!(K, domain::QuasiStaticDomain, Uu::V, U, ka_backend::KernelAbstractionsBackend) where V <: AbstractVector{<:Number}
  update_fields!(U, domain, Uu, ka_backend)
  stiffness!(K, domain, U, ka_backend)
end

function stiffness(domain::QuasiStaticDomain, Uu::V, ka_backend::KernelAbstractionsBackend) where V <: AbstractVector{<:Number}
  U = Cthonios.create_fields(domain, eltype(Uu), ka_backend)
  update_fields!(U, domain, Uu, ka_backend)
  update_bcs!(U, domain.coords, domain.time.current_time, domain.bcs, ka_backend)
  return stiffness(domain, U, ka_backend)
end

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

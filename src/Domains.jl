abstract type Domain <: AbstractCthoniosType end

struct QuasiStaticDomain{
  Coords    <: NodalField,
  Dof       <: DofManager,
  Sections  <: NamedTuple,
  Assembler <: StaticAssembler,
  BCs       <: NamedTuple,
  Time      <: TimeStepper,
  Post      <: PostProcessor,
  #
  Disp      <: NodalField,
  Unknown   <: AbstractArray{<:Number, 1}
} <: Domain
  coords::Coords
  dof::Dof
  sections::Sections
  assembler::Assembler
  bcs::BCs
  time::Time
  post_processor::Post
  #
  U::Disp
  Uu::Unknown
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
  post_processor        = PostProcessor(mesh_file, "output.e", FiniteElementContainers.num_dimensions(mesh_file) |> Int64) # TODO change this
  U                     = FiniteElementContainers.create_fields(dof)
  update_unknown_ids!(dof, bcs)
  Uu                    = FiniteElementContainers.create_unknowns(dof)
  update_bcs!(U, coords, time, bcs)
  FiniteElementContainers.update_fields!(U, dof, Uu)

  # do first assembly pass just to initialize stuff
  for section in sections
    assemble!(assembler, section, coords, U)
  end

  return QuasiStaticDomain(coords, dof, sections, assembler, bcs, time, post_processor, U, Uu)
end

bounbdary_conditions(d::QuasiStaticDomain) = d.bcs
coordinates(d::QuasiStaticDomain)          = d.coords
dof_manager(d::QuasiStaticDomain)          = d.dof
sections(d::QuasiStaticDomain)             = d.sections

num_dimensions(d::QuasiStaticDomain) = FiniteElementContainers.num_fields(coordinates(d))

create_fields(d::QuasiStaticDomain)   = FiniteElementContainers.create_fields(d.dof)
create_fields(d::QuasiStaticDomain, type) = FiniteElementContainers.create_fields(d.dof, type)
create_unknowns(d::QuasiStaticDomain) = FiniteElementContainers.create_unknowns(d.dof)
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

function update_bcs!(d::QuasiStaticDomain)
  for bc in d.bcs
    for node in bc.nodes
      # get coords
      X = @views d.coords[:, node]
      t = d.time.current_time
      d.U[bc.dof, node] = bc.func(X, t)
    end
  end
end

update_fields!(domain::QuasiStaticDomain) = FiniteElementContainers.update_fields!(domain.U, domain.dof, domain.Uu)
update_fields!(U::V1, domain::QuasiStaticDomain, Uu::V2) where {V1 <: NodalField, V2 <: AbstractArray{<:Number, 1}} =
FiniteElementContainers.update_fields!(U, domain.dof, Uu)
update_unknowns!(U::V1, domain::QuasiStaticDomain, Uu::V2) where {V1 <: NodalField, V2 <: AbstractArray{<:Number, 1}} =
FiniteElementContainers.update_fields!(U, domain.dof, Uu)
################################################################

# useful for adjoint calculations and other AD
function energy(domain::QuasiStaticDomain, X::V1, U::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}
  W = 0.0 # Unitful error here # TODO 
  for section in domain.sections
    W = W + energy(section, X, U)
  end
  W
end

function energy(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector{<:Number}
  # Need to update BCs or AD won't do the right thing
  FiniteElementContainers.update_fields!(domain.U, domain.dof, Uu)
  update_bcs!(domain.U, domain.coords, domain.time, domain.bcs)
  return energy(domain, domain.coords, domain.U)
end

function residual_ad end
function energy_func end
function energy_func! end
function energy_gradient_ad end
function energy_hvp_ad end
function hvp end
function jvp end

function energy_gradient(domain::QuasiStaticDomain, X::V, U::V) where V <: NodalField
  domain.assembler.R .= 0.0 # Unitful issue here
  for section in domain.sections
    # residual!(domain.assembler.R, section, X, U)
    residual!(domain.assembler, section, X, U)
  end
  # return domain.assembler.R[domain.dof.is_unknown]
  # return domain.assembler.R

  # R = NodalField{size(X, 1), size(X, 2), Vector, eltype(X)}(undef)
  # R .= domain.assembler.R
  R = NodalField{size(X, 1), size(X, 2), Vector}(domain.assembler.R)
  return R
end

function residual(domain::QuasiStaticDomain, X::V1, U::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}
  domain.assembler.R .= 0.0 # Unitful issue here
  for section in domain.sections
    # residual!(domain.assembler.R, section, X, U)
    residual!(domain.assembler, section, X, U)
  end
  return domain.assembler.R[domain.dof.is_unknown]
end

function residual_for_AD(domain::QuasiStaticDomain, X::V1, U::V2) where {V1 <: AbstractArray, V2 <: AbstractArray}
  # R = Cthonios.create_unknowns(domain, eltype(U)) # Unitful issue here
  R = Cthonios.create_fields(domain, eltype(U))
  for section in domain.sections
    residual!(R, section, X, U)
  end
  return R[domain.dof.is_unknown]
end

function residual(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector{<:Number}
  FiniteElementContainers.update_fields!(domain.U, domain.dof, Uu)
  update_bcs!(domain.U, domain.coords, domain.time, domain.bcs)
  return residual(domain, domain.coords, domain.U)
end

function residual_for_AD(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector
  U = Cthonios.create_fields(domain, eltype(Uu))
  update_fields!(U, domain, Uu)
  update_bcs!(U, domain.coords, domain.time, domain.bcs)
  return residual_for_AD(domain, domain.coords, U)
end

function stiffness(domain::QuasiStaticDomain, X::V, U::V) where V <: NodalField
  domain.assembler.K .= 0.0 # Unitful issue here
  for section in domain.sections
    stiffness!(domain.assembler.K, section, X, U)
  end
  K = domain.assembler.K[domain.dof.is_unknown, domain.dof.is_unknown]
  return 0.5 * (K + K')
end

function stiffness(domain::QuasiStaticDomain, Uu::V) where V <: AbstractVector{<:Number}
  FiniteElementContainers.update_fields!(domain.U, domain.dof, Uu)
  update_bcs!(domain.U, domain.coords, domain.time, domain.bcs)
  return stiffness(domain, domain.coords, domain.U)
end

function assemble!(domain::QuasiStaticDomain, X::V, U::V) where V <: NodalField
  domain.assembler.R .= 0.0
  domain.assembler.K .= 0.0
  for section in domain.sections
    assemble!(domain.assembler, section, X, U)
  end
end

function jvp(domain::QuasiStaticDomain, X::V, U::V, v::A) where {V <: NodalField, A <: AbstractArray{<:Number, 1}}
  # domain.assembler.R .= 0.0
  # domain.assembler.K .= 0.0
  @assert length(v) == length(domain.dof.unknown_indices)
  v_n = NodalField{size(U, 1), size(U, 2), Vector, eltype(U)}(undef)
  v_n .= 0.0 # Unitful issue
  y_n = zeros(Float64, size(U, 1) * size(U, 2))
  FiniteElementContainers.update_fields!(v_n, domain.dof, v)
  for section in domain.sections
    jvp!(y_n, section, X, U, v_n)
  end
  return y_n[domain.dof.is_unknown]
end

# more convenient methods when not worrying about AD or adjoints
energy(domain::QuasiStaticDomain)    = energy(domain, domain.coords, domain.U)
residual(domain::QuasiStaticDomain)  = residual(domain, domain.coords, domain.U)
stiffness(domain::QuasiStaticDomain) = stiffness(domain, domain.coords, domain.U)
assemble!(domain::QuasiStaticDomain) = assemble!(domain, domain.coords, domain.U)

jvp(domain::QuasiStaticDomain, v::V) where V <: AbstractArray{<:Number, 1} = 
jvp(domain, domain.coords, domain.U, v)

##########################################################################

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

function read_bcs(input_settings::D, mesh_file) where D <: Dict
  new_section("Boundary Conditions")
  @info "Reading Boundary Conditions"
  @warn "Currently only supporting displacement bcs"
  @assert "displacement" in keys(input_settings)
  
  @info "Reading Displacement BC nodesets from Exodus file"
  # bcs = EssentialBC[]
  bcs = Dict()
  nset_ids = nodeset_ids(mesh_file)
  ND = FiniteElementContainers.num_dimensions(mesh_file) |> Int64
  for (n, bc) in enumerate(input_settings["displacement"])
    @assert "nodeset id" in keys(bc)
    @assert "dof" in keys(bc)
    @assert bc["nodeset id"] in nset_ids

    nset_id    = bc["nodeset id"]
    nset_name  = read_name(mesh_file.mesh_obj, NodeSet, nset_id)
    if nset_name == ""
      nset_name = "cthonios_nodelist_$nset_id"
    end
    dof        = bc["dof"]
    nset_nodes = convert.(Int64, nodeset(mesh_file, nset_id))
    @show size(nset_nodes)
    nset_nodes = NodalField{1, length(nset_nodes), Vector}(nset_nodes)

    @info "  Displacement BC $n"
    @info "    Nodeset ID   = $nset_id"
    @info "    Nodeset Name = $nset_name"
    @info "    Dof          = $dof"

    if "function" in keys(bc)
      @info "    Function     = $(bc["function"])"
      bc = DisplacementBC{ND, Float64, Float64}(nset_nodes, dof, bc["function"])
    else
      bc = DisplacementBC{ND, Float64, Float64}(nset_nodes, dof)
    end
    bcs[Symbol("displacement_bc_$(nset_name)_dof_$(dof)")] = bc
    @info ""
  end 

  return NamedTuple(bcs)
end

function read_materials(input_settings::D) where D <: Dict
  new_section("Material Models")
  models = Dict{String, ConstitutiveModels.ConstitutiveModel}()
  props  = Dict{String, Vector{Float64}}()
  states = Dict{String, Any}()
  for mat_name in keys(input_settings)
    @assert "model" in keys(input_settings[mat_name])
    @assert "properties" in keys(input_settings[mat_name])

    model_name = input_settings[mat_name]["model"]
    props_in   = input_settings[mat_name]["properties"]

    @info "Reading material $mat_name"
    @info "  Model      = $model_name"
    @info "  Properties = "
    for (key, val) in props_in
      @info "    $(rpad(key, 20)) = $val"
    end
    @info ""

    model, prop, state = eval(Meta.parse(model_name))(props_in)

    models[mat_name] = model
    props[mat_name]  = prop
    states[mat_name] = state
  end
  return models, props, states
end

function read_sections(input_settings, mesh, dof, models, props, states)
  new_section("Sections")
  block_ids = element_block_ids(mesh)
  sections = Dict()
  n = 1
  for section in input_settings
    @info "Reading Section $n"
    @warn "Defaulting to fully integrated element, e.g. q_degree = 2"
    @assert "block id" in keys(section)
    @assert "formulation" in keys(section)
    @assert "material" in keys(section)

    block_id    = section["block id"]
    formulation = section["formulation"]
    mat_name    = section["material"]
    q_degree    = 2 # TODO make this input somehow

    @assert block_id in block_ids
    @assert mat_name in keys(models)
    @assert mat_name in keys(props)
    @assert mat_name in keys(states)

    @info "  Block id      = $block_id"
    @info "  Formulation   = $formulation"
    @info "  Material name = $mat_name"
    @info ""

    conns     = element_connectivity(mesh, block_id)
    conns     = convert.(Int64, conns)
    conns     = Connectivity{size(conns, 1), size(conns, 2), Vector}(conns)
    
    elem_type = FiniteElementContainers.element_type(mesh, block_id)
    if formulation == "plane strain"
      form = FiniteElementContainers.PlaneStrain()
    elseif formulation == "three dimensional" 
      form = FiniteElementContainers.ThreeDimensional()
    else
      @assert false "Unsupported type"
    end

    # fspace  = NonAllocatedFunctionSpace(conns, dof, block_id, q_degree, elem_type)
    fspace  = NonAllocatedFunctionSpace(dof, conns, q_degree, elem_type)
    section = Section(fspace, form, models[mat_name], props[mat_name], states[mat_name]) 

    sections[Symbol("block_$block_id")] = section
    n = n + 1
  end

  # end_section("Sections")

  return NamedTuple(sections)
  # return NamedTuple{Tuple(keys(sections))}(values(sections))
end

###################################################################


struct QuasiStaticDomainOperator{T}
  domain::QuasiStaticDomain
end

Base.eltype(::QuasiStaticDomainOperator{T}) where T = T
Base.size(A::QuasiStaticDomainOperator) = 
(A.domain.dof.unknown_indices |> length, A.domain.dof.unknown_indices |> length)
Base.size(A::QuasiStaticDomainOperator, ::Int) = A.domain.dof.unknown_indices |> length
function Base.:*(A::QuasiStaticDomainOperator, v)
  return jvp(A.domain, v)
end
function LinearAlgebra.mul!(y, A::QuasiStaticDomainOperator, v)
  y = jvp(A.domain, v)
end

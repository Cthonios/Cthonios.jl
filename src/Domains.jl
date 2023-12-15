abstract type AbstractDomain <: AbstractCthoniosType end

struct Domain{
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
}
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

function Domain(input_settings::D) where D <: Dict
  # TODO make variable maybe?
  # arr_type = Vector
    
  for key in required_domain_keys
    @assert key in keys(input_settings)
  end

  # break up input settings
  mesh_settings = input_settings["mesh"]
  bc_settings   = input_settings["boundary conditions"]
  mat_settings  = input_settings["materials"]
  sec_settings  = input_settings["sections"]
  time_settings = input_settings["time stepper"]

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
  Uu                    = FiniteElementContainers.create_unknowns(dof)

  return Domain(coords, dof, sections, assembler, bcs, time, post_processor, U, Uu)
end

bounbdary_conditions(d::Domain) = d.bcs
coordinates(d::Domain)          = d.coords
dof_manager(d::Domain)          = d.dof
sections(d::Domain)             = d.sections

num_dimensions(d::Domain) = FiniteElementContainers.num_fields(coordinates(d))

create_fields(d::Domain)   = FiniteElementContainers.create_fields(d.dof)
create_unknowns(d::Domain) = FiniteElementContainers.create_unknowns(d.dof)

# TODO place where units won't work
function update_unknown_ids!(d::Domain)

  dof = dof_manager(d)
  for bc in bounbdary_conditions(d)
    FiniteElementContainers.update_unknown_ids!(dof_manager(d), bc.nodes, bc.dof)
  end
  
  resize!(dof.unknown_indices, sum(dof.is_unknown))
  dof.unknown_indices .= FiniteElementContainers.dof_ids(dof)[dof.is_unknown]
end

function update_bcs!(d::Domain)
  for bc in d.bcs
    for node in bc.nodes
      # get coords
      X = @views d.coords[:, node]
      t = d.time.current_time
      d.U[bc.dof, node] = bc.func(X, t)
    end
  end
end

function energy(domain::Domain, X::V, U::V) where V <: NodalField
  W = 0.0 # Unitful error here # TODO 
  for section in domain.sections
    W = W + energy(section, X, U)
  end
  W
end

# function energy_grad!(domain::Domain, X::V, dX::V, U::V, dU::V) where V <: NodalField
#   autodiff_deferred(Reverse, energy, Active, domain, Duplicated(X, dX), Duplicated(U, dU))
#   nothing
# end

# function energy_grad_U!(domain::Domain, X::V, U::V, dU::V) where V <: NodalField
#   autodiff_deferred(Reverse, energy, Active, domain, X, Duplicated(U, dU))
#   nothing
# end

# function energy_grad_X!(domain::Domain, X::V, dX::V, U::V) where V <: NodalField
#   autodiff_deferred(Reverse, energy, Active, domain, Duplicated(X, dX), U)
#   nothing
# end

# function residual_ad(domain::Domain, X::V, U::V) where V <: NodalField
#   dU = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   dU .= 0.0
#   energy_grad_U!(domain, X, U, dU)
#   return dU[domain.dof.is_unknown]
# end

function residual(domain::Domain, X::V, U::V) where V <: NodalField
  domain.assembler.R .= 0.0 # Unitful issue here
  for section in domain.sections
    residual!(domain.assembler.R, section, X, U)
  end
  return domain.assembler.R[domain.dof.is_unknown]
end

function stiffness(domain::Domain, X::V, U::V) where V <: NodalField
  domain.assembler.K .= 0.0 # Unitful issue here
  for section in domain.sections
    stiffness!(domain.assembler.K, section, X, U)
  end
  return domain.assembler.K[domain.dof.is_unknown, domain.dof.is_unknown]
end

function solve(domain::Domain)
  reset!(domain.time)

  time_step = 1
  write_time(domain.post_processor, time_step, 0.0)
  write_values(domain.post_processor, NodalVariable, time_step, "displ_x", domain.U[1, :])
  write_values(domain.post_processor, NodalVariable, time_step, "displ_y", domain.U[2, :])

  while domain.time.current_time <= domain.time.end_time
    # time
    step!(domain.time)

    # update bcs
    update_unknown_ids!(domain)
    resize!(domain.Uu, sum(domain.dof.is_unknown))
    update_bcs!(domain)
    update_fields!(domain.U, domain.dof, domain.Uu)

    # some temp calcs
    # W    = energy(domain, domain.coords, domain.U)
    # R_ad = residual_ad(domain, domain.coords, domain.U)

    # TODO wrap in solver
    for n in 1:10
      update_fields!(domain.U, domain.dof, domain.Uu)
      R   = residual(domain, domain.coords, domain.U)
      K   = stiffness(domain, domain.coords, domain.U)
      # @show det(K)

      ΔUu = K \ R
      domain.Uu .= domain.Uu .- ΔUu

      norm_R = norm(R)
      norm_U = norm(ΔUu)

      @info "Iteration $n: ||R|| = $(norm_R)    ||ΔUu|| = $(norm_U)"

      if norm_R < 1e-12 || norm_U < 1e-12
        break
      end
    end

    # post processing
    time_step += 1
    update_fields!(domain.U, domain.dof, domain.Uu)
    write_time(domain.post_processor, time_step, domain.time.current_time)
    write_values(domain.post_processor, NodalVariable, time_step, "displ_x", domain.U[1, :])
    write_values(domain.post_processor, NodalVariable, time_step, "displ_y", domain.U[2, :])
  end

  close(domain.post_processor)
end

# function vp_U!(domain::Domain, X::V1, U::V1, dU::V1, v::V2) where {V1 <: NodalField, V2 <: AbstractVector}
#   @assert length(domain.dof.unknown_indices) == length(v)
#   energy_grad_U!(domain, X, U, dU)
#   R = dU[domain.dof.is_unknown]
#   dot(R, v)
# end

# function residual_ad(domain::Domain, X::V, U::V) where V <: NodalField
#   dU = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   dU .= 0.0
#   energy_grad_U!(domain, X, U, dU)
#   return dU[domain.dof.is_unknown]
# end

# function jvp_U!(domain::Domain, X::V1, U::V1, v::V2) where {V1 <: NodalField, V2 <: AbstractVector{<:Number}}
#   dU = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   d2U = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   d3U = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   dU .= 0.0; d2U .= 0.0; d3U .= 0.0
#   autodiff(Forward, vp_U!, domain, X, Duplicated(U, dU), Duplicated(d2U, d3U), v)
#   return d3U[domain.dof.is_unknown]
# end

# # residual(domain::Domain::U)

# function residual!(R::V, domain::Domain, Uu::V) where V <: AbstractVector{<:Number}
#   update_fields!(domain.U, domain.dof, Uu)
#   R_temp = residual(domain, domain.coords, domain.U)
#   R     .= R_temp[domain.dof.is_unknown]
# end

# function jvp_U!(y::V, domain::Domain, Uu::V, v::V) where V <: AbstractVector{<:Number}
#   update_fields!(domain.U, domain.dof, Uu)
#   y = jvp_U!(domain, domain.coords, domain.U, v)
# end

# function load_step_old!(domain::Domain)
#   # update time
#   step!(domain.time)

#   # update bcs
#   update_unknown_ids!(domain)
#   resize!(domain.Uu, sum(domain.dof.is_unknown))
#   update_bcs!(domain)

#   # make a linear map
#   lm = LinearMap()

# end

# function load_step_old!(domain::Domain)
#   # update time
#   step!(domain.time)

#   # update bcs
#   update_unknown_ids!(domain)
#   resize!(domain.Uu, sum(domain.dof.is_unknown))
#   update_bcs!(domain)

#   # do a solve here to update displacement
#   # dUu = zeros(eltype(domain.Uu), size(domain.Uu))
#   v  = ones(eltype(domain.Uu), size(domain.Uu))
#   dU = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   # d2U = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   # d3U = NodalField{size(domain.U, 1), size(domain.U, 2), Vector, Float64}(undef)
#   # dX = NodalField{size(domain.coords, 1), size(domain.coords, 2), Vector, Float64}(undef)
#   # dU .= 0.0; d2U .= 0.0; d3U .= 0.0
#   # @show dUu
#   W   = energy(domain, domain.coords, domain.U)
#   R_u = energy_grad_U!(domain, domain.coords, domain.U, dU)
#   # vp  = vp_U(domain, domain.coords, domain.U, dU, v)
#   jvp = jvp_U(domain, domain.coords, domain.U, v)
#   R   = residual(domain, domain.coords, domain.U)
#   dU
#   # display(jvp)
#   # R_x = energy_grad_X(domain, domain.coords, dX, domain.U)
#   # R   = energy_grad(domain, domain.coords, dX, domain.U, dU)
#   # dU
#   # display(dU |> sum)
#   # display(dX)
# end

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


# struct DomainOperator{T}
#   domain::Domain
# end

# Base.eltype(::DomainOperator{T}) where T = T
# Base.size(A::DomainOperator) = (A.domain.dof.unknown_indices |> length, A.domain.dof.unknown_indices |> length)
# Base.size(A::DomainOperator, ::Int) = A.domain.dof.unknown_indices |> length
# function Base.:*(A::DomainOperator, v)
#   y = create_unknowns(A.domain)
#   jvp_U!(A.domain, A.domain.coords, A.domain.U, v)
# end
# function LinearAlgebra.mul!(y, A::DomainOperator, v)
#   y = jvp_U!(A.domain, A.domain.coords, A.domain.U, v)
# end


# function solve_domain(domain::Domain)
#   reset!(domain.time)
#   A = DomainOperator{Float64}(domain)

#   # TODO wrap in loop

#   # update time
#   step!(domain.time)

#   # update bcs
#   update_unknown_ids!(domain)
#   resize!(domain.Uu, sum(domain.dof.is_unknown))
#   update_bcs!(domain)

#   # R   = residual(domain, domain.coords, domain.U)
#   # Uu  = domain.Uu
#   # v   = create_unknowns(domain)
#   # y   = create_unknowns(domain)

#   for n in 1:1000
#     # update_fields!(domain.U, domain.dof, domain.Uu)
#     R  = residual(domain, domain.coords, domain.U)
#     domain.Uu .= domain.Uu .- 1e-3 * R
#     update_fields!(domain.U, domain.dof, domain.Uu)

#     @show n norm(R) 
#   end

#   # jvp_U!(y, domain, Uu, v)
  
#   # sol = IterativeSolvers.cg(A, R)
#   # A * v

#   # A * v + y
#   # mul!(y, A, v)
#   # size(A, 1)
#   # eltype(A)
# end
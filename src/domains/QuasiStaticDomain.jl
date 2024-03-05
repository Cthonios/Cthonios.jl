struct QuasiStaticDomainCache{Time, V1, V2, V3, V4, V5, V6, V7} <: AbstractDomainCache
  time::Time
  X::V1
  Uu::V2
  ΔUu::V2
  U::V3
  props::V4
  state_old::V5
  state_new::V5
  Π::V6
  Πs::V7 # for energy kernels
  f::V3  # for internal force
  V::V3 # for working with krylov methods
  Hv::V3 # for working with krylov methods
end

function Base.similar(cache::QuasiStaticDomainCache)
  @unpack time, X, Uu, ΔUu, U, props, state_old, state_new, Π, Πs, f, V, Hv = cache
  return QuasiStaticDomainCache(
    similar(time),
    similar(X), similar(Uu), similar(ΔUu), similar(U), 
    similar(props),
    similar(state_old), similar(state_new),
    similar(Π), similar(Πs), similar(f), similar(V), similar(Hv)
  )
end

function unpack(cache::QuasiStaticDomainCache)
  @unpack time, X, Uu, ΔUu, U, props, state_old, state_new, props, Π, Πs, f, V, Hv = cache
  return time, X, Uu, ΔUu, U, props, state_old, state_new, Π, Πs, f, V, Hv
end

struct QuasiStaticDomain{
  Dof,
  Funcs,
  BCNodes,
  BCDofs,
  BCFuncIDs,
  Sections,
  DomainCache
} <: AbstractDomain{Dof, Funcs, BCNodes, BCDofs, BCFuncIDs, Sections}
  dof::Dof
  funcs::Funcs
  bc_nodes::BCNodes
  bc_dofs::BCDofs
  bc_func_ids::BCFuncIDs
  sections::Sections
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
  sections, props, state_old, state_new, Πs = setup_sections(get_domain_sections_inputs(input_settings), mesh_file, dof)
  # time stepper setup
  time = ConstantTimeStepper(get_domain_time_stepper_inputs(input_settings))

  # cache setup
  Uu = FiniteElementContainers.create_unknowns(dof)
  ΔUu = FiniteElementContainers.create_unknowns(dof)
  U = FiniteElementContainers.create_fields(dof)
  f = FiniteElementContainers.create_fields(dof)
  V = FiniteElementContainers.create_fields(dof)
  Hv = FiniteElementContainers.create_fields(dof)
  domain_cache = QuasiStaticDomainCache(time, coords, Uu, ΔUu, U, props, state_old, state_new, zeros(Float64, 1), Πs, f, V, Hv)

  return QuasiStaticDomain(
    dof, funcs, 
    disp_bc_nodes, disp_bc_dofs, disp_bc_func_ids,
    sections, domain_cache
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
  print(io, "        Type = $(typeof(domain.domain_cache.time).name.name)")
end

# access methods
coordinates(d::QuasiStaticDomain)  = d.coords
dof_manager(d::QuasiStaticDomain)  = d.dof
sections(d::QuasiStaticDomain)     = d.sections
# time_stepper(d::QuasiStaticDomain) = d.time

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
  resize!(d.domain_cache.Uu, length(d.dof.unknown_dofs))
  resize!(d.domain_cache.ΔUu, length(d.dof.unknown_dofs))
end

function update_bcs!(U, domain::QuasiStaticDomain, Xs)
  for (node, dof, func_id) in zip(domain.bc_nodes, domain.bc_dofs, domain.bc_func_ids)
    X = @views Xs[:, node]
    t = domain.domain_cache.time.current_time
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
step!(domain::QuasiStaticDomain) = step!(domain.domain_cache.time)

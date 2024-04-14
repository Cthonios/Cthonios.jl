"""
$(TYPEDFIELDS)
$(TYPEDSIGNATURES)
Section for total Lagrangia (the actual way to spell his name) formulations
"""
struct TotalLagrangeSection{
  ID, FS, Form, Mod
} <: Section{ID, FS, Form, Mod}
  block_id::ID
  fspace::FS
  formulation::Form
  model::Mod
end

function Base.show(io::IO, section::TotalLagrangeSection)
  print(io, "          TotalLagrangeSection\n",
        "            Block id                  = $(section.block_id)\n",
        "            Function space            = $(typeof(section.fspace).name.name)\n",
        "            Formulation               = $(section.formulation)\n",
        "            Material model            = $(section.model)\n")
end

"""
TODO we can remove DofManager as an input if we modify
  the NonAllocatedFunctionSpace constructors
"""
function TotalLagrangeSection(
  mesh::FileMesh,
  dof::DofManager,
  model, 
  formulation,
  block_id::Int,
  q_degree::Int
)
  # setup up conn for fspace
  conns = convert.(Int64, element_connectivity(mesh, block_id))
  conns = Connectivity{size(conns), Vector}(conns)

  # get element type
  elem_type = FiniteElementContainers.element_type(mesh, block_id)

  # setup function space
  fspace = NonAllocatedFunctionSpace(dof, conns, q_degree, elem_type)

  section = TotalLagrangeSection(block_id, fspace, formulation, model)
  return section
end
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

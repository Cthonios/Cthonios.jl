struct TotalLagrangeSection{
  FS, Form, Mod
} <: Section{FS, Form, Mod}
  fspace::FS
  formulation::Form
  model::Mod
end

function TotalLagrangeSection(
  fspace::F, formulation::Form, model::M
) where {
  F <: FunctionSpace, Form <: AbstractMechanicsFormulation,
  M <: MechanicalModel
}

  return TotalLagrangeSection{
    typeof(fspace), typeof(formulation),
    typeof(model)
  }(fspace, formulation, model)
end

function Base.show(io::IO, section::TotalLagrangeSection)
  print(io, "          TotalLagrangeSection\n",
        "            Function space            = $(typeof(section.fspace).name.name)\n",
        "            Formulation               = $(section.formulation)\n",
        "            Material model            = $(section.model)\n")
end

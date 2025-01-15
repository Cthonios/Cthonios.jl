using CairoMakie
using Cthonios
using Documenter
using Literate
using Meshes

DocMeta.setdocmeta!(Cthonios, :DocTestSetup, :(using Cthonios); recursive=true)

function update_mesh_location(content)
  content = replace(content, "Base.source_dir() * \"" => "\"../../../examples/hole_array")
  return content
end

LITERATE_OUTPUT = joinpath(@__DIR__, "src/generated/")

# Literate.markdown(
#   joinpath(@__DIR__, "../examples/hole_array/script.jl"), 
#   LITERATE_OUTPUT,
#   name="hole_array";
#   preprocess=update_mesh_location
# )

makedocs(;
  # modules=[Cthonios],
  authors="Craig M. Hamel <cmhamel32@gmail.com> and contributors",
  repo="https://github.com/Cthonios/Cthonios.jl/blob/{commit}{path}#{line}",
  source="src",
  sitename="Cthonios.jl",
  format=Documenter.HTML(;
    repolink="https://github.com/Cthonios/Cthonios.jl",
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://cthonios.github.io/Cthonios.jl",
    edit_link="main",
    assets=String[],
    size_threshold_warn=5 * 102400,
    size_threshold=10 * 102400
    # size_threshold_ignore=[
    #   "./generated/hole_array.md"
    # ],
  ),
  pages=[
    "Home"                => "index.md",
    "Running Cthonios"    => "running_cthonios.md",
    "Boundary Conditions" => "bcs.md",
    "Domains"             => "domains.md",
    "Iterators"           => "iterators.md",
    "Materials"           => "materials.md",
    "Objectives"          => "objectives.md",
    "Physics"             => "physics.md",
    "Sections"            => "sections.md",
    "Solvers"             => [
      "Linear solvers"    => "linear_solvers.md"
      "Nonlinear solvers" => "nonlinear_solvers.md"
    ],
    "Examples"            => [
      "Hole Array" => "./generated/hole_array.md"
    ]
    # "TimeSteppers"        => "time_steppers.md"
  ],
)

deploydocs(;
  repo="github.com/Cthonios/Cthonios.jl",
  devbranch="main",
)

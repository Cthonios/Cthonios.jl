using Cthonios
using Documenter

DocMeta.setdocmeta!(Cthonios, :DocTestSetup, :(using Cthonios); recursive=true)

makedocs(;
    # modules=[Cthonios],
    authors="Craig M. Hamel <cthonios32@gmail.com> and contributors",
    repo="https://github.com/Cthonios/Cthonios.jl/blob/{commit}{path}#{line}",
    source="src",
    sitename="Cthonios.jl",
    format=Documenter.HTML(;
        repolink="https://github.com/Cthonios/Cthonios.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cthonios.github.io/Cthonios.jl",
        edit_link="main",
        assets=String[],
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
        "TimeSteppers"       => "time_steppers.md"
    ],
)

deploydocs(;
    repo="github.com/Cthonios/Cthonios.jl",
    devbranch="main",
)

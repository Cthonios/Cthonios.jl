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
        "Home"              => "index.md",
        "Running Cthonios"  => "running_cthonios.md",
        "Domains"           => "domains.md",
        "Functions"         => "functions.md",
        "Materials"         => "materials.md",
        "Sections"          => "sections.md",
        "Nonlinear Solvers" => "nonlinear_solvers.md",
        "Linear Solvers"    => "linear_solvers.md"
    ],
)

deploydocs(;
    repo="github.com/Cthonios/Cthonios.jl",
    devbranch="main",
)

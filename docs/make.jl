using Cthonios
using Documenter

DocMeta.setdocmeta!(Cthonios, :DocTestSetup, :(using Cthonios); recursive=true)

makedocs(;
    modules=[Cthonios],
    authors="Craig M. Hamel <cthonios32@gmail.com> and contributors",
    repo="https://github.com/Cthonios/Cthonios.jl/blob/{commit}{path}#{line}",
    sitename="Cthonios.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cthonios.github.io/Cthonios.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Cthonios/Cthonios.jl",
    devbranch="main",
)

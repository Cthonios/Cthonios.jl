using Cthonios
using Documenter

DocMeta.setdocmeta!(Cthonios, :DocTestSetup, :(using Cthonios); recursive=true)

makedocs(;
    modules=[Cthonios],
    authors="Craig M. Hamel <cmhamel32@gmail.com> and contributors",
    repo="https://github.com/cmhamel/Cthonios.jl/blob/{commit}{path}#{line}",
    sitename="Cthonios.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cmhamel.github.io/Cthonios.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cmhamel/Cthonios.jl",
    devbranch="main",
)

using MusicSpiralRepresentation
using Documenter

DocMeta.setdocmeta!(MusicSpiralRepresentation, :DocTestSetup, :(using MusicSpiralRepresentation); recursive=true)

makedocs(;
    modules=[MusicSpiralRepresentation],
    authors="'Alfredo González-Espinoza'",
    repo="https://github.com/spiralizing/MusicSpiralRepresentation.jl/blob/{commit}{path}#{line}",
    sitename="MusicSpiralRepresentation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://spiralizing.github.io/MusicSpiralRepresentation.jl",
        edit_link="main",
        assets=String[],
        mathengine=Documenter.MathJax3(),
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
    ],
)

deploydocs(;
    repo="github.com/spiralizing/MusicSpiralRepresentation.jl",
    devbranch="main",
)

using SphericalGraphEmbedding
using Documenter

DocMeta.setdocmeta!(SphericalGraphEmbedding, :DocTestSetup, :(using SphericalGraphEmbedding); recursive=true)

makedocs(;
    modules=[SphericalGraphEmbedding],
    authors="Michaël Fanuel <mrfanuel@hotmail.fr> and contributors",
    repo="https://github.com/mrfanuel/SphericalGraphEmbedding.jl/blob/{commit}{path}#{line}",
    sitename="SphericalGraphEmbedding.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mrfanuel.github.io/SphericalGraphEmbedding.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mrfanuel/SphericalGraphEmbedding.jl",
    devbranch="main",
)

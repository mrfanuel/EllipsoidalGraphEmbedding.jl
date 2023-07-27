using EllipsoidalGraphEmbedding
using Documenter

DocMeta.setdocmeta!(EllipsoidalGraphEmbedding, :DocTestSetup, :(using EllipsoidalGraphEmbedding); recursive=true)

makedocs(;
    modules=[EllipsoidalGraphEmbedding],
    authors="Michaël Fanuel <mrfanuel@hotmail.fr> and contributors",
    repo="https://github.com/mrfanuel/EllipsoidalGraphEmbedding.jl/blob/{commit}{path}#{line}",
    sitename="EllipsoidalGraphEmbedding.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mrfanuel.github.io/EllipsoidalGraphEmbedding.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mrfanuel/EllipsoidalGraphEmbedding.jl",
    devbranch="main",
)

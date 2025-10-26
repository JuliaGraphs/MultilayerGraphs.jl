using Revise
using SimpleTraits
using Distributions
using Graphs
using MultilayerGraphs
using Documenter

DocMeta.setdocmeta!(
    MultilayerGraphs, :DocTestSetup, :(using MultilayerGraphs); recursive=true
)

makedocs(;
    modules=[MultilayerGraphs],
    authors="Pietro Monticone, Claudio Moroni",
    sitename="MultilayerGraphs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliagraphs.org/MultilayerGraphs.jl",
        assets=String[],
        size_threshold_ignore=["API.md"],
    ),
    pages=["🏠 Home" => "index.md", "🛠 API" => "API.md"],
    clean=false,
    warnonly=true,
)

deploydocs(;
    repo="github.com/JuliaGraphs/MultilayerGraphs.jl", devbranch="main", push_preview=true
)

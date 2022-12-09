using Revise
using MultilayerGraphs
using Documenter

# cd("dev/MultilayerGraphs/docs")

DocMeta.setdocmeta!(
    MultilayerGraphs, :DocTestSetup, :(using MultilayerGraphs); recursive=true
)

makedocs(;
    modules=[MultilayerGraphs],
    authors="Pietro Monticone, Claudio Moroni",
    repo="https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/{commit}{path}#{line}",
    sitename="MultilayerGraphs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliagraphs.org/MultilayerGraphs.jl",
        assets=String[],
    ),
    pages=["Home" => "index.md", "API" => "API.md"],
    clean = false
)

deploydocs(;
    repo="github.com/JuliaGraphs/MultilayerGraphs.jl", devbranch="main", push_preview=true
)

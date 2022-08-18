using MultilayerGraphs
using Documenter

DocMeta.setdocmeta!(
    MultilayerGraphs, :DocTestSetup, :(using MultilayerGraphs); recursive=true
)

makedocs(;
    modules=[MultilayerGraphs],
    authors="Pietro Monticone, Claudio Moroni",
    repo="https://github.com/InPhyT/MultilayerGraphs.jl/blob/{commit}{path}#{line}",
    sitename="MultilayerGraphs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://InPhyT.github.io/MultilayerGraphs.jl",
        assets=String[],
    ),
    pages=["Home" => "index.md", "Internals" => "internals.md"],
)

deploydocs(;
    repo="github.com/InPhyT/MultilayerGraphs.jl", devbranch="main", push_preview=true
)

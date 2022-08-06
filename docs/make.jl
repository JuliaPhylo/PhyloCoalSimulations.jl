using PhyloCoalSimulations
using Documenter

DocMeta.setdocmeta!(PhyloCoalSimulations, :DocTestSetup, :(using PhyloCoalSimulations); recursive=true)

makedocs(;
    modules=[PhyloCoalSimulations],
    authors="Cecile Ane <cecileane@users.noreply.github.com> and contributors",
    repo="https://github.com/cecileane/PhyloCoalSimulations.jl/blob/{commit}{path}#{line}",
    sitename="PhyloCoalSimulations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cecileane.github.io/PhyloCoalSimulations.jl",
        assets=String[],
    ),
    pages=[
        "home" => "index.md",
        "manual" => [
            "getting started" => "man/getting_started.md",
            "more examples" => "man/mapping_genetree_to_network.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/cecileane/PhyloCoalSimulations.jl",
    devbranch="main",
    push_preview = true,
)

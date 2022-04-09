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
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cecileane/PhyloCoalSimulations.jl",
    devbranch="main",
)

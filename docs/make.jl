using PhyloCoalSimulations
using Documenter

DocMeta.setdocmeta!(PhyloCoalSimulations, :DocTestSetup, :(using PhyloCoalSimulations); recursive=true)

makedocs(;
    modules=[PhyloCoalSimulations],
    authors="Cecile Ane <cecileane@users.noreply.github.com> and contributors",
    sitename="PhyloCoalSimulations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cecileane.github.io/PhyloCoalSimulations.jl/stable/",
        assets=String[],
    ),
    pages=[
        "home" => "index.md",
        "manual" => [
            "getting started" => "man/getting_started.md",
            "mapping gene trees into the species network" => "man/mapping_genetree_to_network.md",
            "converting between units" => "man/converting_coal2generation_units.md",
            "correlated inheritance" => "man/correlated_inheritance.md",
            "more examples" => "man/more_examples.md",
        ],
        "library" => [
            "public" => "lib/public.md",
            "internals" => "lib/internal.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/cecileane/PhyloCoalSimulations.jl",
    devbranch="main",
    push_preview = true,
)

# notes to maintain documentation

- built with [Documenter](https://juliadocs.github.io/Documenter.jl)
- deployed [here](https://cecileane.github.io/PhyloCoalSimulations.jl/)
  (go to `dev/` or `stable/`)
  using github and files committed to the `gh-pages` branch.

The workflow is controlled by `.github/workflows/documentation.yml`,
which launches `docs/make.jl`.
To meet the dependency on R (via PhyloPlots) documentation.yml installs R.

## what to update

- Documenter version in `docs/Project.toml`
- update Julia version in `.github/workflows/documentation.yml`

## to build the documentation locally

Run the commands below to:
- test the `jldoctest` blocks of examples in the docstrings
- create or update a `build/` directory with html files

```shell
julia --project=docs/ -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs/ --color=yes docs/make.jl
```

or interactively in `docs/`:
```shell
pkg> activate .
pkg> instantiate
pkg> dev ~/.julia/dev/PhyloCoalSimulations
julia> include("make.jl")
```

Then open `src/build/index.html`.

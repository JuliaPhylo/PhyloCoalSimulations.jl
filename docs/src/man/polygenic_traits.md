```@setup polygenictrait
using PhyloNetworks, PhyloCoalSimulations, StableRNGs
```

# polygenic traits

We can now simulate traits controlled by loci that evolve according to
the coalescent.

[`simulate_polygenictrait`](@ref) simulates polygenic traits,
controlled by multiple loci, whose effects on the trait are additive,
that is:
``X = (Y_1+...+ Y_L)/\sqrt{L}``
for trait X is controlled by L loci.
Each locus effect Y evolves under its own gene tree.
The same process is used across all loci.
For each locus, each ancestor in the root population has its effect
drawn from some "root" distribution. The locus effect evolves according
to some "transition" distribution along each edge of its gene tree.

The documentation for [`simulate_polygenictrait`](@ref) shows two
examples, one with a normal distribution at the root of the species
phylogeny and normal transition probabilities (like a Brownian motion)
for each locus on its gene tree.
The other example shows 0/1 binary traits (e.g. major/minor allele)
with a Markov transition probability.

The example below shows how to build a compound Poisson process,
to model a mutational model for each locus. Under this model,
the effect of each locus stays constant until a mutation affects
the locus. Mutations occur with a Poisson process at some rate.
At each mutation, the locus effect changes by some value
(mutational effect) drawn from some distribution.

We build such a process below, from which we can draw samples
along one edge of a gene tree.

```@example polygenictrait
using Distributions, Random
import Distributions: rand

"""
    CompoundPoissonGaussian(m, μ, σ)

Compound Poisson distribution with Gaussian-distributed effects.
`m` is the mean, or the starting value for the process.
`μ` is the average number of mutations, sampled from a Poisson distribution.
`σ` is the standard deviation of each mutation effect, sampled from a
  Gaussian distribution (of mean 0, standard deviation `σ`).

Example: we may define `d = CompoundPoissonGaussian(3, 0.6, 1)`
for the evolution of a locus effect starting at value 3,
evolving along an edge whose length predicts 0.6 mutations,
and with a standard deviation of 1 for each mutation effect.
Then we can sample the locus effect value at the end of the edge as:
`rand(d)`.
We could also draw 4 samples, say, with `rand(d, 4)`.
"""
struct CompoundPoissonGaussian <: Sampleable{Univariate,Continuous}
    "mean"
    m::Float64
    "mean number of events, sampled from Poisson(μ)"
    μ::Float64
    "standard deviation of each event effect, sampled from Normal(0,σ)"
    σ::Float64
end
function rand(rng::AbstractRNG, s::CompoundPoissonGaussian)
    n = rand(rng, Poisson(s.μ)) # number of mutations
    return s.m + sum(rand(rng, Normal(0,s.σ), n))
end
nothing # hide
```

We can now use this transition distribution to simulate
a polygenic trait under a mutational process.
We use the same network as used earlier.

```@example polygenictrait
net = readnewick("((C:0.9,(B:0.2)#H1:0.7::0.6)i1:0.6,(#H1:0.6::0.4,A:1.0)i2:0.5)i3;");
nothing # hide
```
![3-taxon network](../assets/figures/net3taxa.svg)

Next we pick parameters:
`μ=1.5` mutations per coalescent unit (aka, per 2N generations where
2N is the haploid effective population size),
`σ=2` for the mutation effects,
and the equilibrium variance at the root population: `μσ²`

```@repl polygenictrait
μ=1.5; σ=2; rootdist = Normal(0, sqrt(μ)*σ);
transition(y0, len) = CompoundPoissonGaussian(y0, μ*len, σ);
```

finally, we simulate 2 replicates of the process, using L=4 loci,
and 2 individuals in each population:

```@repl polygenictrait
x,lab,gt = simulate_polygenictrait(net, 2, 4, rootdist, transition;
    nindividuals=2);
lab
```

The trait values in `x` come in the same order as in the labels `lab`.

```@repl polygenictrait
using DataFrames
pop = replace.(lab, r"_\d+$" => "")
DataFrame(individual=lab, population=pop, rep1=x[1][:], rep2=x[2][:])
```

In case it's useful, the gene trees are returned,
with the L=4 gene trees from the first replicated listed first,
the L trees used for the second replicate listed next, etc.
```@repl polygenictrait
length(gt) # 8 gene trees: 4 for replicate 1, then 4 for replicate 2
writenewick(gt[1], round=true)
```

If one does *not* wish to save the gene trees, we could tell julia
like this below (with the dash),
so that memory holding gene trees can be recycled sooner:

```julia
x,lab,_ = simulate_polygenictrait(...);
```

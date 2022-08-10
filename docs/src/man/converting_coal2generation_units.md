```@setup converting
using PhyloNetworks, PhyloCoalSimulations
net = readTopology("((C:0.9,(B:0.2)#H1:0.7::0.6)I1:0.6,(#H1:0.6::0.4,A:1.0)I2:0.5)I3;");
using Random; Random.seed!(261); # as in mapping block
tree = simulatecoalescent(net,1,1; nodemapping=true)[1];
```
# converting coalescent units to number of generations

Edge lengths in gene trees are simulated in coalescent units.
These lengths can be converted into numbers of generations by multiplying by
the effective population size Nₑ, since coalescent units are `u = g/Nₑ`.
This can be done with different Nₑ's across different edges in the network,
including the root edge.

!!! info "diploid versus haploid Nₑ"
    The formula `u = g/Nₑ` uses the haploid effective population size Nₑ.
    For diploid taxa and autosomes, the haploid population size should be twice
    the diploid population size, for example.

Here is an example using the same network and simulated gene tree as earlier.
```@repl converting
writeTopology(net)
writeTopology(tree, round=true)
```
![example 1, same as in mapping section](../assets/figures/genetree_example1.svg)

Let's set Nₑ set to 1,000 in all populations (including the population above the root),
except in edge 6. For this edge 6 (population leading to species A),
let's set its Nₑ to 10,000.

```@repl converting
Ne = Dict(e.number => 1_000 for e in net.edge);
push!(Ne, 8 => 1_000); # add Ne for the edge above the network's root
Ne[6] = 10_000;        # higher population size for the edge to species A
Ne
writeTopology(tree, round=true) # lengths in coalescent units: before unit conversion
# convert edge lengths in gene tree from coalescent units to # generations
for e in tree.edge
  e.length = round(e.length * Ne[e.inCycle]) # round: to get integers
end
writeTopology(tree, round=true) # lengths in # of generations
```

Note that the simulation model assumes an infinite-Nₑ approximation,
so the rescaling of edge lengths from coalescent units to number of generations
will be imperfect for very small populations size. With the extreme Nₑ=1,
coalescences should be immediate in a single generation back in time: g=1.
Using the approximation, the simulated number of generations will typically be
between 0-3 generations. But this is an extreme case, and the approximation
should be very good even for moderate Nₑ's.


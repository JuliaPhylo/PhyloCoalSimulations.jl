```@setup getting_started
using PhyloNetworks, PhyloCoalSimulations, RCall
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
using Random; Random.seed!(432)
```

# getting started

## installation

To install Julia see [here](https://docs.julialang.org/en/v1/manual/getting-started/)
and to install Julia packages, see [here](https://docs.julialang.org/en/v1/stdlib/Pkg/).
To install `PhyloCoalSimulations` in the Julia REPL
(as well as `PhyloNetworks` for many utilities), enter package mode with `]`
and do this:

```
add PhyloCoalSimulations
add PhyloNetworks
```
or do this in julian mode:

```julia
using Pkg
Pkg.add("PhyloCoalSimulations")
Pkg.add("PhyloNetworks")
```

## basic simulation example

### example network

For a basic example, we use a simple 3-species network plotted below.
On the left, the plot shows the edge numbers (black) and the γ inheritance values (blue).
On the right, the length of horizontal lines are proportional to edge lengths,
and the plot shows the edge length values.

```@example getting_started
using PhyloNetworks
net = readTopology("((C:0.9,(B:0.2)#H1:0.7::0.6):0.6,(#H1:0.6,A:1):0.5);");
using PhyloPlots
R"svg"(figname("net3taxa.svg"), width=6, height=3); # hide
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2]); # hide
plot(net, :R, showEdgeNumber=true, showIntNodeLabel=true, showGamma=true, tipOffset=0.1);
R"mtext"("in black: edge numbers", side=1, line=-1);  # hide
plot(net, :R, showEdgeLength=true, useEdgeLength=true, tipOffset=0.1);
R"mtext"("in black: edge lengths", side=1, line=-1);  # hide
R"dev.off()" # hide
nothing # hide
```
![3-taxon network](../assets/figures/net3taxa.svg)

Note that this example network is not time consistent: the length of the path
from the root to the hybridization node H1 is different depending if we go
through the major edge (0.6+0.7=1.3) or the minor edge (0.5+0.6=1.1).

Coalescent simulations can be performed along such networks, also
along non-ultrametric networks.
If the network is ultrametric (time-consistent, and with all tips at the
same distance from the root), then gene trees will also be ultrametric.

### basic example: simulate, save to file, plot

We use [`simulatecoalescent`](@ref) to simulate gene trees along this network.
Below, we simulate 2 gene trees. By default, there's 1 individual per species.

```@repl getting_started
trees = simulatecoalescent(net, 2)
```

Branch lengths are assumed to be in coalescent units in the species network
(number of generations / effective population size), and edge lengths in gene
trees are also in coalescent units.

We can work with these gene trees within Julia with downstream code,
and/or we can save them to a file:

```@repl getting_started
writeMultiTopology(trees, stdout) # write them to standout output (screen here)
writeMultiTopology(trees, "genetrees.phy") # warning: will overwrite "genetrees.phy" if this file existed
rm("genetrees.phy") # hide
```

Let's plot these 2 gene trees. In the plots below, we annotate each
edge with its attribute that tracked the network edge on which
the coalescent event occured (where the gene tree lineage originated,
going back in time).
For example, the gene lineage that ends in A is always mapped to network edge
6, which is the number of the external edge to A in the network (see plot
of network [above](#example-network) on the left).

```@example getting_started
R"svg"(figname("genetrees_gettingstarted_1.svg"), width=6, height=3); # hide
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2]); # hide
using DataFrames
for i in 1:2
  gt = trees[i]
  plot(gt, :R, tipOffset=0.1, useEdgeLength=true,
               edgeLabel=DataFrame(number = [e.number  for e in gt.edge],
                                   label  = [e.inCycle for e in gt.edge]));
  R"mtext"("gene $i", line=-1) # hide
end
R"mtext"("numbers: network edge each gene lineage maps to, at time of coalescence.\n8 = number of edge above the network root", side=1, line=-1, outer=true);  # hide
R"dev.off()" # hide
nothing # hide
```
![3-taxon network](../assets/figures/genetrees_gettingstarted_1.svg)

### several individuals per species

We can ask for more individuals. To simulate 3 individuals / species for example:

```@repl getting_started
simulatecoalescent(net, 1, 3) # 1 gene tree only. 3 individuals in each species
```

We can also ask for varying numbers of individuals. For example, we simulate
below 2 individuals in A and 1 individual in each of B and C,
using a dictionary to map species to their number of individuals:

```@repl getting_started
genetrees = simulatecoalescent(net, 1, Dict("A"=>2, "B"=>1, "C"=>1));
writeTopology(genetrees[1])
```

We can set 0 individuals within a species to simulate missing data.

```@repl getting_started
genetrees = simulatecoalescent(net, 3, Dict("A"=>2, "B"=>1, "C"=>0));
writeMultiTopology(genetrees, stdout)
```

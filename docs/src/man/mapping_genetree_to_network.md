```@setup mapping
using PhyloNetworks, PhyloPlots, PhyloCoalSimulations, RCall, DataFrames
figname(x) = joinpath("..", "assets", "figures", x)
```
# mapping gene trees into the species phylogeny

Nodes in gene trees can be mapped to a node or an edge in the species phylogeny.
Edges in gene trees can be mapped to an edge in the species phylogeny.
Attributes of nodes and edges are used to carry this mapping information,
as detailed in the documentation of function [`simulatecoalescent`](@ref).

We give examples below of how we may use this mapping information.

## naming internal nodes in the network and gene trees

First, it's useful to name internal nodes in the network, to which we
can later map nodes in the gene tree.

```@repl mapping
net = readTopology("((C:0.9,(B:0.2)#H1:0.7::0.6):0.6,(#H1:0.6,A:1):0.5);");
PhyloNetworks.nameinternalnodes!(net, "I"); # "I" is a prefix to name internal nodes
writeTopology(net)
```

Notice the extra node names in the species phylogeny: I1, I2, and I3 at the root.

Next, we use the option `nodemapping=true` when simulating gene trees,
to ask for extra degree-2 nodes in gene trees. These nodes are created each time
that a gene tree lineage crosses a node in the species phylogeny, that is,
each time that a gene tree lineage crosses a speciation node or a hybridization node.
We'll simulate a single tree here each time.

```@repl mapping
using Random; Random.seed!(261); # to replicate randomness
tree_regular = simulatecoalescent(net,1,1)[1];
writeTopology(tree_regular, round=true) # regular nodes only
Random.seed!(261); # to replicate the same coalescent simulation
tree = simulatecoalescent(net,1,1; nodemapping=true)[1];
writeTopology(tree, round=true) # extra degree-2 nodes for mapping
```

Notice that regular nodes in the gene tree (nodes that we get without the
`nodemapping` option) don't have names. With the `nodemapping` option, there are
many new nodes, all of degree-2, named after internal nodes in the network
(I1, I2, I3). These extra degree-2 nodes and their names are sufficient to map
the gene tree into the species network.

The network is shown on the left below, with edges annotated by their numbers.

```@example mapping
R"svg"(figname("genetrees_moreexample1.svg"), width=6, height=3); # hide
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2]); # hide
plot(net, :R, showEdgeNumber=true, showIntNodeLabel=true, tipOffset=0.1);
R"mtext"("species network", side=3, line=-1);  # hide
plot(tree, :R, edgeLabel=DataFrame(number=[e.number for e in tree.edge],
                                   label=[e.inCycle for e in tree.edge]),
               showIntNodeLabel=true, tipOffset=0.1);
R"mtext"("gene tree", side=3, line=-1);  # hide
R"dev.off()" # hide
nothing # hide
```
![example 1: degree-2 node names in gene tree](../assets/figures/genetrees_moreexample1.svg)

In the gene tree (right), each lineage is annotated by the network
edge it maps into. Degree-2 nodes appear via their names, such that each
horizontal line represents a series of gene lineages, separated from each other
by degree-2 nodes.
For example, the horizontal line tracing B's ancestry back in time maps into the
network like this:
- from B, go back along edge 2
- meet hybrid node H1, was inherited from the minor hybrid edge 5,
- from speciation node I2, trace back along edge 7,
- meet the network's root node I3, trace back along the network's root edge 8
  before coalescing with the ancestor of the other lineages (which have already
  coalesced by then).

## cleaning gene trees

Almost all examples below use this mapping information via the extra degree-2
nodes and the extra edges between these nodes.

But we may want to "clean" gene trees of their degree-2 nodes at some point.
This can be done with the `PhyloNetworks` utility `removedegree2nodes!`, like this:

```julia
PhyloNetworks.removedegree2nodes!(tree, true)
```
The option `true` is to keep the root, even if it's of degree 2.

!!! note
    PhyloNetworks v0.15.0 does not have this option and
    removes the root of degree 2, but the option will be available in the next
    version of PhyloNetworks.
    This is not a problem for simulating sequences along gene trees when using a
    reversible substitution model, for which the root placement doesn't matter.

# converting coalescent units to number of generations

Edge lengths in gene trees are simulated in coalescent units.
These lengths can be converted into numbers of generations by multiplying by
the effective population size Nₑ, since coalescent units are `u = g/Nₑ`.
This can be done with different Nₑ's across different edges in the network,
including the root edge. Here is an example, with Nₑ set to 1,000 in all
populations (including the population above the root),
except in edge 6: the population leading to species A has its Nₑ set to 10,000.

```@repl mapping
Ne = Dict(e.number => 1_000 for e in net.edge);
push!(Ne, 8 => 1_000); # add Ne for the edge above the network's root
Ne[6] = 10_000;        # higher population size for the edge to species A
Ne
# convert edge lengths in gene tree from coalescent units to # generations
tree_in_generations = deepcopy(tree)
for e in tree_in_generations.edge
  e.length = round(e.length * Ne[e.inCycle]) # round: to get integers
end
writeTopology(tree_in_generations, round=true, digits=4) # after rate variation
```

Note that the simulation model assumes an infinite-Nₑ approximation,
so the rescaling of edge lengths from coalescent units to number of generations
will be imperfect for very small populations size. With the extreme Nₑ=1,
coalescences should be immediate in a single generation back in time: g=1.
Using the approximation, the simulated number of generations will typically be
between 0-3 generations. But this is an extreme case, and the approximation
should be very good even for moderate Nₑ's.

# example uses

##  counting deep coalescences

The number of deep coalescences can be quantified as the number of
"extra" lineages due to incomplete lineage sorting, that can be calculated
from embedding the gene tree into the species phylogeny
(see [Maddison 1997](https://doi.org/10.1093/sysbio/46.3.523) for species trees).
For a speciation node in the species tree (e.g. I1), there are extra lineages
in the gene tree, owing to a lack of coalescence, if there are more than 2 gene
lineages mapping to this speciation node.
For each hybridization node (e.g. H1), there are extra lineages if there are
more than 1 gene lineage mapping to this hybridization node
(because there's only 1 child edge descending from a hybridization node).

We can count the number of extra lineages by counting the number of degree-2
nodes in the gene tree mapping to each node in the species network, then
counting how many are "extra".
In our gene tree above, we can do it this way:

```@repl mapping
node_count = Dict(n.name => 0      for n in net.node if !n.leaf) # ignore leaves
node_hyb = Dict(n.name => n.hybrid for n in net.node if !n.leaf) # true/1 for hybrid nodes
# traverse the gene tree, to count number of lineages entering each network node
for n in tree.node
  (n.leaf || n.name == "") && continue # skip leaves and nodes without a name
  node_count[n.name] += 1  # increment by 1 the number of lineages entering this node
end
node_count # the 2 lineages entering I2 didn't coalesce until they reached I3
for node in keys(node_count) # modify our counts, to only keep the extras
  node_count[node] = max(0, node_count[node] - ( node_hyb[node] ? 1 : 2))
end
node_count # number of extra lineages: 1 extra entering I3
deepcoalescence = sum(values(node_count))
```
On the particular gene tree we simulated, we counted 1 deep coalescence.

## number of lineages inherited via gene flow

Our network has inheritance γ=0.4 on the minor edge, which we'll call the
"gene flow" edge, and γ=0.6 on the major hybrid edge, parent to H1 on the major tree.
But we may be interested in the realized proportion of lineages inherited
from each parent at H1, realized in the gene trees we actually simulated.
To do so, we can count the number of gene lineages that are mapped to each
hybrid edge in the network.

This mapping is stored in the edge attribute `.inCycle`.
From the plot above, the minor "gene flow" edge is edge number 5 and the
major hybrid edge has number 3.
So we can count the gene lineages inherited via gene flow
as the number of gene tree edges with `inCycle` equal to 5.

If the gene trees have been saved to a file and later read from this file,
then the `.inCycle` attributes are no longer stored in memory. In this case,
we can retrieve the mapping information by the internal node names.
The edges going through gene flow are those whose child node is named "H1"
and parent node is named "I2".

We use the first option with the `.inCycle` attribute below.
We get that our one simulated gene tree was indeed inherited via gene flow:

```@repl mapping
sum(e.inCycle == 5 for e in tree.edge)
```

To make this more interesting, we can simulate 100 gene trees
then count how many were inherited via gene flow. If we ask for 2 individuals
in species B, then each gene may have 2 lineages that enter the hybrid node H1,
if the two B individuals fail to coalesce. In that case, it's possible that one
individual lineage was inherited via gene flow, and the other not.
We'll calculate the gene flow proportion among all these lineages.
This proportion should be close (but not exactly equal) to the theoretical
γ=0.4 from the network model.

```@repl mapping
ngenes = 100
genetrees = simulatecoalescent(net, ngenes, Dict("B"=>2, "A"=>1, "C"=>1); nodemapping=true);
length(genetrees)
nlineages_geneflow = sum(sum(e.inCycle == 5 for e in gt.edge) for gt in genetrees)
nlineages_major    = sum(sum(e.inCycle == 3 for e in gt.edge) for gt in genetrees)
proportion_geneflow = nlineages_geneflow / (nlineages_geneflow + nlineages_major)
```

## rate variation across species

The gene trees resulting from `simulatecoalescent` have their edge lengths
in coalescent units. One may want to convert them to substitutions per site,
so as to simulate molecular sequences along these gene trees.
The mapping information is important to allow for different rates of molecular
evolution across different species. Here is an example to do this.

We will use the `Distributions` package to simulate rates from a log-normal
distribution across species, that is, across edges in the species network.

```@repl mapping
using Distributions
lognormal_rate_dist = LogNormal(-0.125, 0.5) # μ = -σ²/2 to get a mean of 1.
networkedge_rate = Dict(e.number => rand(lognormal_rate_dist) for e in net.edge)
# add entry for the edge above the network's root. Find its number first.
rootedgenumber = maximum(e.number for e in net.edge) + 1
push!(networkedge_rate, rootedgenumber => rand(lognormal_rate_dist))
writeTopology(tree, round=true, digits=4) # before rate variation
# multiply the length of each gene lineage by the rate of the species edge it maps into
for e in tree.edge
  e.length *= networkedge_rate[e.inCycle]
end
writeTopology(tree, round=true, digits=4) # after rate variation
```

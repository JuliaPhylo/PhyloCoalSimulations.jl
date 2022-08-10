```@setup downstreamexamples
using PhyloNetworks, PhyloCoalSimulations
net = readTopology("((C:0.9,(B:0.2)#H1:0.7::0.6)I1:0.6,(#H1:0.6::0.4,A:1.0)I2:0.5)I3;");
using Random; Random.seed!(261); # as in mapping block
tree = simulatecoalescent(net,1,1; nodemapping=true)[1];
```
# example uses

We are re-using the same network and simulated tree as before:
```@repl downstreamexamples
writeTopology(net)
writeTopology(tree, round=true)
```
![example 1, same as in mapping section](../assets/figures/genetree_example1.svg)

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

```@repl downstreamexamples
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

```@repl downstreamexamples
sum(e.inCycle == 5 for e in tree.edge)
```

To make this more interesting, we can simulate many gene trees
then count how many of their lineages were inherited via gene flow.
If we ask for 2 individuals in species B,
then each gene may have 2 lineages that enter the hybrid node H1,
if the two B individuals fail to coalesce. In that case, it's possible that one
individual lineage was inherited via gene flow, and the other not.
We'll calculate the gene flow proportion among all these lineages.
This proportion should be close (but not exactly equal) to the theoretical
γ=0.4 from the network model.

```@repl downstreamexamples
ngenes = 100
genetrees = simulatecoalescent(net, ngenes, Dict("B"=>2, "A"=>1, "C"=>1); nodemapping=true);
length(genetrees)
nlineages_geneflow = sum(sum(e.inCycle == 5 for e in gt.edge) for gt in genetrees)
nlineages_major    = sum(sum(e.inCycle == 3 for e in gt.edge) for gt in genetrees)
# realized γ, close to 0.4:
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

```@repl downstreamexamples
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

```@setup downstreamexamples
using PhyloNetworks, PhyloCoalSimulations, StableRNGs
net = readnewick("((C:0.9,(B:0.2)#H1:0.7::0.6)i1:0.6,(#H1:0.6::0.4,A:1.0)i2:0.5)i3;");
using Random; Random.seed!(261) # for examples below that use default RNG
rng = StableRNG(7); # as in mapping block
tree = simulatecoalescent(rng, net,1,1; nodemapping=true)[1];
```
# example uses

We are re-using the same network and simulated tree as before:
```@repl downstreamexamples
writenewick(net)
writenewick(tree, round=true)
```
![example 1, same as in mapping section](../assets/figures/genetree_example1.svg)

##  counting deep coalescences

The number of deep coalescences can be quantified as the number of
"extra" lineages due to incomplete lineage sorting, that can be calculated
from embedding the gene tree into the species phylogeny
(see [Maddison 1997](https://doi.org/10.1093/sysbio/46.3.523) for species trees).
For an edge in the network, say edge 7 going from i2 to i3 going in back in time,
lineage sorting is complete if all the gene lineages entering the edge (at i2)
coalesce into a single gene lineage by the time they exit the edge (at i3).
If they don't, the number of extra lineages is `k-1` where `k` is the number of
lineages "exiting" the edge, for that particular edge in the species network
and that particular gene tree.
The total number of extra lineages, for a given gene tree, is the sum across
all edges in the species phylogeny.

In our gene tree above, we can count the number of lineages that exit each
species edge using the degree-2 mapping nodes,
then count how many lineages are "extra". We do so below using utilities
[`mappingnodes`](@ref PhyloCoalSimulations.mappingnodes)
to iterate over degree-2 mapping nodes and
[`population_mappedto`](@ref)
to extract the mapping information.

```@repl downstreamexamples
# dictionary to store the count of extra lineages exiting each network edge. initialized at 0
edge_count = Dict(e.number => 0 for e in net.edge)
const PCS = PhyloCoalSimulations; # for lazy typing!
for n in PCS.mappingnodes(tree)  # iterate over degree-2 mapping nodes in the gene tree
  child = getchildedge(n)
  popid = population_mappedto(child) # number of species edge that 'n' came from
  # sanity check below
  isnothing(popid) && error("""population ID not found for the child edge of
                    node number $(n.number) mapping to species node $(n.name).""")
  edge_count[popid] += 1 # increment by 1 the number of lineages exiting population edge 'popid'
end
edge_count
```

From this, we see two interesting things.
- 0 lineages exited edge number 3 in the species network: it's the hybrid
  edge from H1 to i3 (going back in time). That's because the only lineage
  at H1 was interited from i2, so there weren't any lineage evolving through edge 3.
- 2 lineages exited edge number 7 (going from i2 to i3 back in time),
  so that's 1 extra lineage. All other edges look as expected, with a single
  gene lineage exiting from them.

We can now calculate the total number of extra lineages:

```@repl downstreamexamples
filter!(p -> p.second > 0, edge_count) # filter out edges without any gene lineage
map!(k -> k-1, values(edge_count))     # calculate number of "extras": k-1
deepcoalescence = sum(values(edge_count)) # sum 'extras' over all edges in the network
```

On the particular gene tree we simulated, we counted 1 deep coalescence.

## number of lineages inherited via gene flow

Our network has inheritance γ=0.4 on the minor edge, which we'll call the
"gene flow" edge, and γ=0.6 on the major hybrid edge, parent to H1 on the major tree.
But we may be interested in the realized proportion of lineages inherited
from each parent at H1, realized in the gene trees we actually simulated.
To do so, we can count the number of gene lineages that are mapped to each
hybrid edge in the network.

This mapping is stored in the edge attribute `.inte1` internally,
but it's best to access it via the function [`population_mappedto`](@ref)
(as the internal representation may change).
From the plot above, the minor "gene flow" edge is edge number 5 and the
major hybrid edge has number 3.
So we can count the gene lineages inherited via gene flow
as the number of gene tree edges with `population_mappedto(edge)` equal to 5.

If the gene trees have been saved to a file and later read from this file,
then they no longer have the internal attributes that store the population
each edge arose from.
But we can retrieve this mapping information from the internal node names.
For example, the edges going through gene flow are those whose child node is
named "H1" and parent node is named "i2".
To recover the internal edge attributes, use [`gene_edgemapping!`](@ref) as
described in section
"[mapping gene tree edges back into species edges](@ref)".

We assume below that the internal edge attribute are up-to-date (as when gene
trees were simulated just before).
We get that our one simulated gene tree was indeed inherited via gene flow:

```@repl downstreamexamples
sum(population_mappedto(e) == 5 for e in tree.edge)
```

Next we define a function to do this for any edge, so we can re-use later:
```@repl downstreamexamples
nlineages_through(edgeID, gt) = sum(population_mappedto(e) == edgeID for e in gt.edge);
nlineages_through(5, tree) # same as before: now done via our new function
nlineages_through(3, tree) # lineages that went through edge 3, the major edge.
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
ngenes = 100;
genetrees = simulatecoalescent(net, ngenes, Dict("B"=>2, "A"=>1, "C"=>1); nodemapping=true);
length(genetrees)
nlineages_geneflow = sum(nlineages_through(5,gt) for gt in genetrees)
nlineages_major    = sum(nlineages_through(3,gt) for gt in genetrees)
proportion_geneflow = nlineages_geneflow / (nlineages_geneflow + nlineages_major)
# realized γ, close to 0.4
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
```

Next we want to simulate a rate for the ancestral edge above the network's root.
We find its number first, then add a entry to our dictionary of rates
```@repl downstreamexamples
rootedgenumber = PhyloCoalSimulations.get_rootedgenumber(net)
push!(networkedge_rate, rootedgenumber => rand(lognormal_rate_dist))
writenewick(tree, round=true, digits=4) # before rate variation
```

Finally, we multiply the length of each gene lineage by the rate of
the species edge it maps into:
```@repl downstreamexamples
for e in tree.edge
  e.length *= networkedge_rate[population_mappedto(e)]
end
writenewick(tree, round=true, digits=4) # after rate variation
```

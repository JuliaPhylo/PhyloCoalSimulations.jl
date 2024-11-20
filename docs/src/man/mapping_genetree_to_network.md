```@setup mapping
using PhyloNetworks, PhyloPlots, PhyloCoalSimulations, RCall, DataFrames, StableRNGs
figpath = joinpath("..", "assets", "figures"); mkpath(figpath)
figname(x) = joinpath(figpath, x)
```
# mapping gene trees into the species phylogeny

Nodes in gene trees can be mapped to a node or an edge in the species phylogeny.
Edges in gene trees can be mapped to an edge in the species phylogeny.
Attributes of nodes and edges are used to carry this mapping information,
as detailed in the documentation of function [`simulatecoalescent`](@ref).

We give examples below of how we may use this mapping information.

## naming internal nodes

First, it's useful to name internal nodes in the network, to which we
can later map nodes in the gene tree.

```@repl mapping
net = readnewick("((C:0.9,(B:0.2)#H1:0.7::0.6):0.6,(#H1:0.6,A:1):0.5);");
nameinternalnodes!(net); # default prefix "i" to name internal nodes
writenewick(net)
```

Notice the extra node names in the species phylogeny: i1, i2, and i3 at the root.

Next, we use the option `nodemapping=true` when simulating gene trees,
to ask for extra degree-2 nodes in gene trees. These nodes are created each time
that a gene tree lineage crosses a node in the species phylogeny, that is,
each time that a gene tree lineage crosses a speciation node or a hybridization node.
We'll simulate a single tree here each time.

```@repl mapping
using StableRNGs; rng = StableRNG(7); # to replicate randomness, but default RNG is better
tree_regular = simulatecoalescent(rng, net,1,1)[1];
writenewick(tree_regular, round=true) # regular nodes only
rng = StableRNG(7); # to replicate the same coalescent simulation
tree = simulatecoalescent(rng, net,1,1; nodemapping=true)[1];
writenewick(tree, round=true) # extra degree-2 nodes for mapping
```

Notice that regular nodes in the gene tree (nodes that we get without the
`nodemapping` option) don't have names. With the `nodemapping` option, there are
many new nodes, all of degree-2, named after internal nodes in the network
(i1, i2, i3). These extra degree-2 nodes and their names are sufficient to map
the gene tree into the species network.

The network is shown on the left below, with edges annotated by their numbers.

```@example mapping
R"svg"(figname("genetree_example1.svg"), width=7.5, height=3); # hide
R"par"(mar=[.1,.2,.1,.2], oma=[0,0,0,0.8]); R"layout"([1 2]); # hide
plot(net, showedgenumber=true, shownodelabel=true, tipoffset=0.1);
R"mtext"("species network", side=3, line=-1);  # hide
R"mtext"("grey: population edge number", side=1, line=-1, cex=0.9);  # hide
plot(tree, edgelabel=DataFrame(number=[e.number for e in tree.edge],
                               label=[e.inte1 for e in tree.edge]),
           edgelabelcolor="red4", shownodelabel=true, tipoffset=0.1);
R"mtext"("gene tree", side=3, line=-1);  # hide
R"mtext"("red (edge inte1 value): population a gene edge maps into", side=1, line=-2, cex=0.9); # hide
R"mtext"("black (node names): speciation/reticulation a node maps to", side=1, line=-1, cex=0.9); # hide
R"dev.off()" # hide
nothing # hide
```
![example 1: degree-2 node names in gene tree](../assets/figures/genetree_example1.svg)

In the gene tree (right), each lineage is annotated by the network
edge it maps into. Degree-2 nodes appear via their names, such that each
horizontal line represents a series of gene lineages, separated from each other
by degree-2 nodes.
For example, the horizontal line tracing B's ancestry back in time maps into the
network like this:
- from B, go back along edge 2
- meet hybrid node H1, was inherited from the minor hybrid edge 5,
- from speciation node i2, trace back along edge 7,
- meet the network's root node i3, trace back along the network's root edge 8
  before coalescing with the ancestor of the other lineages (which have already
  coalesced by then).

## cleaning gene trees

Almost all examples below use this mapping information via the extra degree-2
nodes and the extra edges between these nodes.

But we may want to "clean" gene trees of their degree-2 nodes at some point.
This can be done with the `PhyloNetworks` utility `removedegree2nodes!`, like this:

```@repl mapping
PhyloNetworks.removedegree2nodes!(tree, true)
```
The option `true` is to keep the root, even if it's of degree 2.

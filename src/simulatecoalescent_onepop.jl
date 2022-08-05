"""
    simulatecoal_onepopulation!(lineagelist, population_length, nextlineageID,
                               populationID=-1)

Simulate the coalescent process within a single population of length given
in coalescent units, starting from lineages in `lineagelist`. This list should
be a vector of incomplete edges, that is, edges incident to a single node only.

Vector of incomplete edges, whose lengths have been increased, is modified in place.
New nodes and their parent edges are created by coalescent events, numbered
with consecutive integers starting at `nextlineageID`.

Output: nextlineageID, incremented by number of new lineages.

In lineages, edge lengths are also considered in coalescent units.

The newly created nodes and edges have their `.inCycle` attribute set to
`populationID`, so as to track the mapping of gene lineages to populations
in the species phylogeny.

# examples

```jldoctest
julia> PhyloCoalSimulations.simulatecoal_onepopulation!([], 2.0, 1)
1

julia> e1 = PhyloCoalSimulations.initializetip("s","1",1,"",0.1);

julia> e2 = PhyloCoalSimulations.initializetip("s","2",2,"",0.2);

julia> forest = [e1,e2]
2-element Vector{PhyloNetworks.Edge}:
 PhyloNetworks.Edge:
 number:1
 length:0.1
 attached to 1 node(s) (parent first): 1

 PhyloNetworks.Edge:
 number:2
 length:0.2
 attached to 1 node(s) (parent first): 2

julia> using Random; Random.seed!(7690);

julia> PhyloCoalSimulations.simulatecoal_onepopulation!(forest, Inf, 3);

julia> PhyloCoalSimulations.convert2tree!(forest[1].node[1])
PhyloNetworks.HybridNetwork, Rooted Network
2 edges
3 nodes: 2 tips, 0 hybrid nodes, 1 internal tree nodes.
tip labels: s2, s1
(s2:0.502,s1:0.402);

```
"""
function simulatecoal_onepopulation!(lineagelist::AbstractVector,
            poplen::AbstractFloat, nextid::Int, populationid = -1)
    isempty(lineagelist) && return(nextid)
    poplen > 0.0 || return(nextid)
    # at this point: the list has 1 or more lineages
    nlineage = length(lineagelist)
    poplen < Inf || nlineage > 1 ||
        error("there should be 2 or more lineages at the start of the infinite root population")
    timeleft = poplen
    while true
        if nlineage == 1 # then no one can coalesce
            lineagelist[1].length += timeleft
            return nextid # breaks the while loop
        end
        # at this point: 2 or more lineages
        coaltime = rand(Exponential(1/binomial(nlineage,2)))
        toadd = min(coaltime, timeleft)
        for e in lineagelist
            e.length += toadd
        end
        timeleft -= coaltime
        if timeleft <= 0 # break the loop: coalescence goes too far
            return nextid
        end
        # pick at random 2 lineages to merge, into a new node numbered nextid
        drop_index = sample(1:nlineage)
        edge1 = popat!(lineagelist, drop_index) # popat! requires julia v1.5
        drop_index = sample(1:(nlineage-1))
        edge2 = lineagelist[drop_index]
        parentedge = coalescence_edge(edge1,edge2,nextid,populationid)
        lineagelist[drop_index] = parentedge
        nextid += 1
        nlineage -= 1
    end
end

"""
    coalescence_edge(edge1, edge2, number, populationid)

Create a coalescence between edges 1 and 2: with a new parent node
`n` numbered `number` and a new parent edge `e`
above the parent node, of length 0 and numbered `number`.
Both `n.inCycle` and `e.inCycle` are set to `populationid`.
"""
function coalescence_edge(e1,e2,number,populationid)
    parentnode = PN.Node(number, false) # false because tree node (not hybrid)
    parentnode.inCycle = populationid
    push!(e1.node, parentnode) # isChild1 true by default
    push!(e2.node, parentnode)
    push!(parentnode.edge, e1)
    push!(parentnode.edge, e2)

    parentedge = PN.Edge(number, 0.0)   # length 0.0
    parentedge.inCycle = populationid
    push!(parentedge.node, parentnode)
    push!(parentnode.edge, parentedge)
    return parentedge
end

"""
    initializetip(species::AbstractString, individual::AbstractString,
                  number::Integer, delim=""::AbstractString, len=0.0)

Create a leaf node and a pendant edge of length `len`, incident to each other,
both numbered `number`. Return the pendant edge.
The leaf name is made by concatenating `species`, `delim` and `individual`.
"""
function initializetip(species::AbstractString, individual::AbstractString,
                       number::Integer, delim=""::AbstractString,
                       len=0.0::AbstractFloat)
    tipnode = PN.Node(number,false)
    tipnode.leaf = true
    tipnode.name = species * delim * individual
    tipedge = PN.Edge(number, len)
    push!(tipedge.node, tipnode)
    push!(tipnode.edge, tipedge)
    return tipedge
end

"""
    initializetip(speciesnode::Node, nindividuals::Integer,
                  number::Integer, delim, len=0.0)

Vector of pendant leaf edges, with leaves named after `speciesnode`,
and numbered with consecutive number IDs starting at `number`.
If nindividuals is 1, then the leaf name is simply the species name.
Otherwise, then the leaf names include the individual number and
the default delimiter is `_`. For example, if the species name is `s`
then leaf names are: `s_1`, `s_2`, etc. by default.
"""
function initializetip(speciesnode::PN.Node, nindividuals::Integer,
                       number::Integer, delim=nothing,
                       len=0.0::AbstractFloat)
    sname = speciesnode.name
    sname != "" || error("empty name: initializing at a non-leaf node?")
    forest = PN.Edge[]
    if isnothing(delim)
        delim = (nindividuals == 1 ? "" : "_")
    end
    iname(x) = (nindividuals == 1 ? "" : string(x))
    for i in 1:nindividuals
        push!(forest, initializetip(sname, iname(i), number, delim, len))
        number += 1
    end
    return forest
end

"""
    map2population!(forest, population_node, populationid, nextlineageID)

Extend each incomplete edge in the forest with a new degree-2 node `n` and
a new incomplete edge `e`, with the following information to map `n` and `e`
into the species phylogeny:
- `e.inCycle` is set to `populationid`, and
- `n.name` is set to `population_node.name` if this name is non-empty, or
  `string(population_node.number)` otherwise (with any negative sign replaced
  by the string "minus").
`e.number` and `n.number` are set to `nextlineageID`, which is incremented by 1
for each incomplete edge in the forest.

The forest is updated to contain the newly-created incomplete edges,
replacing the old incomplete (and now complete) edges.

Output: nextlineageID, incremented by the number of newly created degree-2 lineages.

# example

```jldoctest
julia> using PhyloNetworks; net = readTopology("(A:1,B:1);");

julia> leafA = net.node[1]; edge2A_number = net.edge[1].number;

julia> f = PhyloCoalSimulations.initializetip(leafA, 2, 4); # 2 edges, numbered 4 & 5

julia> PhyloCoalSimulations.map2population!(f, leafA, edge2A_number, 6)
8

julia> length(f)
2

julia> f[2]
PhyloNetworks.Edge:
 number:7
 length:0.0
 attached to 1 node(s) (parent first): 7

julia> [e.node[1].name for e in f]
2-element Vector{String}:
 "A"
 "A"
```
"""
function map2population!(forest, pop_node, populationid, number)
    popnodename = (pop_node.name == "" ? replace(string(pop_node.number), "-" => "minus") :
                                         pop_node.name)
    for i in eachindex(forest)
        e_old = forest[i]
        degree2node = PN.Node(number, false) # false because tree node (not hybrid)
        # keep default inCycle = -1, because maps to population node, not population edge
        degree2node.name = popnodename
        push!(e_old.node, degree2node) # isChild1 true by default
        push!(degree2node.edge, e_old)
        e_new = PN.Edge(number, 0.0)   # length 0.0
        e_new.inCycle = populationid
        push!(e_new.node, degree2node)
        push!(degree2node.edge, e_new)
        forest[i] = e_new
        number += 1
    end
    return number
end

function cleanrootnode!(rootnode::PN.Node)
    for j in length(rootnode.edge):-1:1
        ee = rootnode.edge[j]
        if length(ee.node) == 1 # incomplete edge: prune it
            deleteat!(ee.node, 1)
            deleteat!(rootnode.edge, j)
        end
    end
end

"""
    convert2tree!(rootnode)

Return a network with all nodes and edges that can be reached from `rootnode`.
**Warning**: Assumes that edges are correctly directed (with correct `isChild1`
attribute) and that the graph is a tree. This is *not* checked.

If the root node is still attached to an incomplete root edge, this
edge & node are first disconnected.
"""
function convert2tree!(rootnode::PN.Node)
    cleanrootnode!(rootnode)
    net = PN.HybridNetwork()
    push!(net.node, rootnode)
    net.root = 1
    net.isRooted = true
    for edge in rootnode.edge
        collect_edges_nodes!(net, edge)
    end
    net.numNodes = length(net.node)
    net.numEdges = length(net.edge)
    length(Set(n.number for n in net.node)) == net.numNodes ||
        @error("node numbers are not all distinct")
    length(Set(e.number for e in net.edge)) == net.numEdges ||
        @error("edge numbers are not all distinct")
    ntaxa = 0
    for nn in net.node # collect leaf names and # of leaves
        nn.leaf == (length(nn.edge)==1) ||
            error("incorrect .leaf for node number $nn")
        nn.leaf || continue
        ntaxa += 1
        nn.name != "" || error("leaf without name")
        push!(net.leaf,  nn)
        push!(net.names, nn.name)
    end
    net.numTaxa = ntaxa
    return net
end
function collect_edges_nodes!(net, parentedge)
    push!(net.edge, parentedge)
    nn = PN.getChild(parentedge)
    push!(net.node, nn)
    for childedge in nn.edge
        childedge !== parentedge || continue # skip parentedge
        collect_edges_nodes!(net, childedge)
    end
    return nothing
end

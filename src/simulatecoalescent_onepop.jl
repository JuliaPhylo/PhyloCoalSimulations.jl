"""
    simulatecoal_onepopulation!(lineagelist, population_length, nextlineageID,
                               populationID=-1)

Simulate the coalescent process within a single population of length given
in coalescent units, starting from lineages in `lineagelist`. This list should
be a vector of incomplete edges, that is, edges incident to a single node only.

Output: vector of incomplete edges, whose lengths have been increased.
New nodes and their parent edges are created by coalescent events, numbered
with consecutive integers starting at `nextlineageID`.

In lineages, edge lengths are also considered in coalescent units.

The newly created nodes and edges have their `.inCycle` attribute set to
`populationID`, so as to track the mapping of gene lineages to populations
in the species phylogeny.

# examples

```jldoctest
julia> PhyloCoalSimulations.simulatecoal_onepopulation!([], 2.0, 1)
Any[]

julia> e1 = PhyloCoalSimulations.initializetip("s","1",1,"",0.1);

julia> e2 = PhyloCoalSimulations.initializetip("s","2",2,"",0.2);

julia> using Random; Random.seed!(7690);

julia> forest = PhyloCoalSimulations.simulatecoal_onepopulation!([e1,e2], Inf, 3);

julia> PhyloCoalSimulations.convert2tree!(forest[1].node[1])
PhyloNetworks.HybridNetwork, Rooted Network
2 edges
3 nodes: 2 tips, 0 hybrid nodes, 1 internal tree nodes.
tip labels: s1, s2
(s2:0.573,s1:0.473);

```
"""
function simulatecoal_onepopulation!(lineagelist::AbstractVector,
            poplen::AbstractFloat, nextid::Int, populationid = -1)
    isempty(lineagelist) && return(lineagelist)
    poplen > 0.0 || return(lineagelist)
    # at this point: the list has 1 or more lineages
    nlineage = length(lineagelist)
    poplen < Inf || nlineage > 1 ||
        error("there should be 2 or more lineages at the start of the infinite root population")
    timeleft = poplen
    while true
        if nlineage == 1 # then no one can coalesce
            lineagelist[1].length += timeleft
            return lineagelist # breaks the while loop
        end
        # at this point: 2 or more lineages
        coaltime = rand(Exponential(1/binomial(nlineage,2)))
        toadd = min(coaltime, timeleft)
        for e in lineagelist
            e.length += toadd
        end
        timeleft -= coaltime
        if timeleft <= 0 # break the loop: coalescence goes too far
            return lineagelist
        end
        # pick at random 2 lineages to merge, into a new node numbered nextid
        drop_index = sample(1:nlineage)
        edge1 = popat!(lineagelist, drop_index)
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
    trackmapping2population(forest, population_nodename, population_edgenumber)

Update each incomplete edge `e` in the forest, to:
- set to `population_edgenumber` the edge's `isCycle` attribute, and
- set to `population_nodename` the name of the first node incident to the edge.

# example

```jldoctest
julia> e1 = PhyloCoalSimulations.initializetip("s","1",1,"",0.1);

julia> PhyloCoalSimulations.trackmapping2population!([e1], "10", 11)

julia> e1.inCycle
11

julia> e1.node[1]
PhyloNetworks.Node:
 number:1
 name:10
 leaf node
 attached to 1 edges, numbered: 1

```
"""
function trackmapping2population!(forest, pop_nodename, pop_edgenumber)
    for e in forest
        e.inCycle = pop_edgenumber
        e.node[1].name = pop_nodename
    end
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
        error("node numbers are not all distinct")
    length(Set(e.number for e in net.edge)) == net.numEdges ||
        error("edge numbers are not all distinct")
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

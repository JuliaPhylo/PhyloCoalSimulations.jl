"""
    simulatecoalescent(net, nloci, nindividuals=1; nodemapping=false)

Simulate `nloci` gene trees with `nindividuals` from each species
under the multispecies network coalescent, along network `net`.
Branch lengths in `net` are interpreted as being in coalescent units
(ratio: number of generations / effective population size).
The coalescent model uses the infinite-population-size approximation.

Output: vector of gene trees, of length `nloci`.

`nindividuals` can be a single integer, or a dictionary listing the
number of individuals to be simulated for each species.

If `nodemapping` is true, each simulated gene tree is augmented with
degree-2 nodes that can be mapped to speciation or hybridization events.
The mapping of gene tree nodes & edges to network edges is carried
by the `.inCycle` attribute. The mapping of gene tree nodes to network nodes
is carried by the `.name` attribute. Namely:

- A degree-3 node (1 parent + 2 children) represents a coalescent event
  that occurred along a population *edge* in `net`. Its `.inCycle` attribute
  is set to the number of that network population edge.
  Its parent edge has its `.inCycle` attribute also set to the number of
  the population edge that it originated from.
- The gene tree's root node (of degree 2) represents a coalescent event
  along the network's root edge.
  Its `.inCycle` attribute is the number assigned to the network's root edge,
  which is set to the maximum edge number + 1.
- A leaf (or degree-1 node) represents an individual. It maps to a species
  in `net`. The individual leaf name is set to the species name
  if `nindividuals` is 1. Otherwise, its name is set to `speciesname_i`
  where `i` is the individual number in that species.
  Its `inCycle` attribute is the default `-1`.
- A non-root degree-2 node represents a speciation or hybridization and maps
  to a population *node* in `net`. Its `inCycle` attribute is the default `-1`.
  Its name is set to network node name, if it exists. If the network node
  has no name, the gene tree node is given a name built from the network node
  number.

Assumptions:
- `net` must have non-missing edge lengths and γ values.
- If `nindividuals` is a dictionary, it must have a key for all species, with
  the same spelling of species names in its keys as in the tip labels of `net`.
# examples

```jldoctest
julia> using PhyloNetworks

julia> net = readTopology("(A:1,B:1);"); # branch lengths of 1 coalescent unit

julia> using Random; Random.seed!(54321); # for replicability of examples below

julia> simulatecoalescent(net, 2) # 2 gene trees, default of 1 individual/species
2-element Vector{HybridNetwork}:
 PhyloNetworks.HybridNetwork, Rooted Network
2 edges
3 nodes: 2 tips, 0 hybrid nodes, 1 internal tree nodes.
tip labels: B, A
(B:1.023,A:1.023);

 PhyloNetworks.HybridNetwork, Rooted Network
2 edges
3 nodes: 2 tips, 0 hybrid nodes, 1 internal tree nodes.
tip labels: B, A
(B:2.328,A:2.328);


julia> simulatecoalescent(net, 1, 3)[1] # 1 gene tree, 3 individuals/species
PhyloNetworks.HybridNetwork, Rooted Network
10 edges
11 nodes: 6 tips, 0 hybrid nodes, 5 internal tree nodes.
tip labels: B_2, B_1, B_3, A_3, ...
(((B_2:0.12,B_1:0.12):2.39,B_3:2.51):0.692,(A_3:0.518,(A_2:0.461,A_1:0.461):0.057):2.684);

julia> simulatecoalescent(net, 1, Dict("A"=>2, "B"=>1))[1] # 2 individuals in A, 1 in B
PhyloNetworks.HybridNetwork, Rooted Network
4 edges
5 nodes: 3 tips, 0 hybrid nodes, 2 internal tree nodes.
tip labels: B, A_2, A_1
(B:2.801,(A_2:0.344,A_1:0.344):2.457);


julia> Random.seed!(21);

julia> tree1 = simulatecoalescent(net,2,2; nodemapping=true)[1]; # first gene tree only

julia> writeTopology(tree1, round=true)
"((B_1:1.0)minus2:0.301,(((A_2:0.11,A_1:0.11):0.89)minus2:0.148,(B_2:1.0)minus2:0.148):0.153);"

julia> PhyloNetworks.nameinternalnodes!(net, "I"); writeTopology(net)
"(A:1.0,B:1.0)I1;"

julia> tree1 = simulatecoalescent(net,2,2; nodemapping=true)[1]; writeTopology(tree1, round=true)
"(((B_2:0.214,B_1:0.214):0.786)I1:0.168,((A_2:1.0)I1:0.008,(A_1:1.0)I1:0.008):0.16);"

julia> printNodes(net)
node leaf  hybrid hasHybEdge name inCycle edges'numbers
1    true  false  false      A    -1      1   
2    true  false  false      B    -1      2   
-2   false false  false      I1   -1      1    2   

julia> printNodes(tree1)
node leaf  hybrid hasHybEdge name inCycle edges'numbers
10   false false  false           3       8    9   
8    false false  false      I1   -1      5    8   
5    false false  false           2       4    3    5   
4    true  false  false      B_2  -1      4   
3    true  false  false      B_1  -1      3   
9    false false  false           3       7    6    9   
7    false false  false      I1   -1      2    7   
2    true  false  false      A_2  -1      2   
6    false false  false      I1   -1      1    6   
1    true  false  false      A_1  -1      1   

julia> [(tree_edge_number = e.number, pop_edge_number = e.inCycle) for e in tree1.edge]
9-element Vector{NamedTuple{(:tree_edge_number, :pop_edge_number), Tuple{Int64, Int64}}}:
 (tree_edge_number = 8, pop_edge_number = 3)
 (tree_edge_number = 5, pop_edge_number = 2)
 (tree_edge_number = 4, pop_edge_number = 2)
 (tree_edge_number = 3, pop_edge_number = 2)
 (tree_edge_number = 9, pop_edge_number = 3)
 (tree_edge_number = 7, pop_edge_number = 3)
 (tree_edge_number = 2, pop_edge_number = 1)
 (tree_edge_number = 6, pop_edge_number = 3)
 (tree_edge_number = 1, pop_edge_number = 1)
```
"""
function simulatecoalescent(net::PN.HybridNetwork, nloci, nindividuals=1;
        nodemapping=false)
    if isa(nindividuals, AbstractDict)
        issubset(PN.tipLabels(net), keys(nindividuals)) || error("nindividuals is missing some species")
        valtype(nindividuals) <: Integer || error("nindividuals should be integers")
    elseif isa(nindividuals, Integer)
        nindividuals = Dict(n.name => nindividuals for n in net.leaf)
    else
        error("nindividuals should be an integer or vector of integers")
    end
    for e in net.edge
        e.hybrid && e.gamma == -1.0 && error("the network needs gamma values")
    end
    PN.preorder!(net)
    node2forest = Dict(n.number => PN.Edge[] for n in net.node)
    edge2forest = Dict(e.number => PN.Edge[] for e in net.edge)
    parentedges = Vector{PN.Edge}[] # parentedge[i] = vector of edges parent to node i in "nodes_changed"
    childedges = Vector{PN.Edge}[]
    nnodes = length(net.node)
    for nodei in 1:nnodes
        nn = net.nodes_changed[nodei]
        parentedgelist = PN.Edge[]
        childedgelist = PN.Edge[]
        for e in nn.edge
            PN.getChild(e) === nn ? push!(parentedgelist, e) : push!(childedgelist, e)
        end
        push!(parentedges, parentedgelist)
        push!(childedges, childedgelist)
    end
    rootedgenumber = max(0, maximum(e.number for e in net.edge)) + 1

    genetreelist = Vector{PN.HybridNetwork}(undef,nloci)
    for ilocus in 1:nloci
        for f in values(node2forest) # clean up intermediate dictionary
            empty!(f)
        end
        for f in values(edge2forest)
            empty!(f)
        end
        nextid = 1 # number for the next node & its parent edge, to be created in gene tree
        for nodei in nnodes:-1:1
            nn = net.nodes_changed[nodei]
            parentedgelist = parentedges[nodei] # parent population(s)
            nparents = length(parentedgelist)
            if nn.leaf
                ee = parentedgelist[1]
                f = initializetip(nn, nindividuals[nn.name], nextid)
                # john: could the inCycle be set inside initializetip?
                for e in f
                    e.inCycle = ee.number
                end
                nextid += nindividuals[nn.name]
                edge2forest[ee.number] = f
                nextid = simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
                continue 
            end
            # gather forest from children
            f = node2forest[nn.number]
            for ce in childedges[nodei] # children edges in network
                append!(f, edge2forest[ce.number])
            end

            if nparents == 1 # tree node but not leaf
                ee = parentedgelist[1]
                if nodemapping # create degree-2 nodes named after nn
                    nextid = map2population!(f, nn, ee.number, nextid)
                end
                edge2forest[ee.number] = f
                nextid = simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
            elseif nparents > 1
                # partition forest's lineages randomly across all parent edges
                d = Distributions.Categorical([e.gamma for e in parentedgelist])
                for genetree in f
                    roll=rand(d)
                    whereto = parentedgelist[roll] # edge in species network
                    push!(edge2forest[whereto.number], genetree)
                end
                # run coalescent along each parent edge
                for ee in parentedgelist
                    f = edge2forest[ee.number]
                    if nodemapping # degree-2 nodes named after nn, and forest edges mapping to chosen ee
                        nextid = map2population!(f, nn, ee.number, nextid)
                    end
                    nextid = simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
                end
            else # nparents = 0: infinite root population
                # f can have a single edge if a displayed tree's MRCA is strictly below the root,
                #   or if the network root has a single child
                if length(f) > 1
                    if nodemapping
                        nextid = map2population!(f, nn, rootedgenumber, nextid)
                    end
                    nextid = simulatecoal_onepopulation!(f, Inf, nextid, rootedgenumber)
                end
                rootnode = f[1].node[1]
                genetreelist[ilocus] = convert2tree!(rootnode)
            end
        end
    end
    return genetreelist
end

"""
    get_rootedgenumber(network)

1 + maximum of 0 and of all the network's edge numbers:
can be used as a unique identifier of the edge above the network's root.
"""
get_rootedgenumber(net) = max(0, maximum(e.number for e in net.edge)) + 1

"""
    simulatecoalescent(net, nloci, nindividuals; nodemapping=false, inheritancecorrelation=0.0)

Simulate `nloci` gene trees with `nindividuals` from each species
under the multispecies network coalescent, along network `net`
whose branch lengths are assumed to be in **coalescent units**
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
  which is set by [`get_rootedgenumber`](@ref) as the maximum edge number + 1.
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

By default, lineages at a hybrid node come from a parent (chosen according
to inheritance probabilities γ) *independently* across lineages.
Positive dependence can be simulated with option `inheritancecorrelation`.
For example, if this correlation is set to 1, then all lineages inherit from
the same (randomly sampled) parent. More generally, the lineages' parents
are simulated according to a Dirichlet process with base distribution determined
by the γ values, and with concentration parameter α = (1-r)/r, that is, r = 1/(1+α),
where `r` is the input inheritance correlation.

Assumptions:
- `net` must have non-missing edge lengths and γ values.
- If `nindividuals` is a dictionary, it must have a key for all species, with
  the same spelling of species names in its keys as in the tip labels of `net`.

# examples

```jldoctest
julia> using PhyloNetworks

julia> net = readTopology("(A:1,B:1);"); # branch lengths of 1 coalescent unit

julia> using Random; Random.seed!(54321); # for replicability of examples below

julia> simulatecoalescent(net, 2, 1) # 2 gene trees, 1 individual/species
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
function simulatecoalescent(net::PN.HybridNetwork, nloci::Integer, nindividuals;
        nodemapping=false, inheritancecorrelation=0.0)
    if isa(nindividuals, AbstractDict)
        issubset(PN.tipLabels(net), keys(nindividuals)) || error("nindividuals is missing some species")
        valtype(nindividuals) <: Integer || error("nindividuals should be integers")
    elseif isa(nindividuals, Integer)
        nindividuals = Dict(n.name => nindividuals for n in net.leaf)
    else
        error("nindividuals should be an integer or dictionary of integers")
    end
    inheritancecorrelation >= 0.0 || error("the inheritance correlation should be non-negative")
    inheritancecorrelation <= 1.0 || error("the inheritance correlation should be <= 1")
    independentlin = (inheritancecorrelation == 0.0)
    alpha = (1-inheritancecorrelation)/inheritancecorrelation # correlation = 1/(1+alpha)
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
    rootedgenumber = get_rootedgenumber(net)

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
                f = initializetipforest(nn, nindividuals[nn.name], nextid)
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

            isempty(f) && continue # to next node

            if nparents == 1 # tree node but not leaf
                ee = parentedgelist[1]
                if nodemapping # create degree-2 nodes named after nn
                    nextid = map2population!(f, nn, ee.number, nextid)
                end
                edge2forest[ee.number] = f
                nextid = simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
            elseif nparents > 1
                # partition forest's lineages randomly across all parent edges
                parentprob = [e.gamma for e in parentedgelist]
                parentdist = Distributions.Categorical(parentprob)
                if independentlin
                    for genetree in f
                      roll = rand(parentdist)
                      push!(edge2forest[parentedgelist[roll].number], genetree)
                    end
                else # then use a Dirichlet process to simulate correlated parents
                    # do one lineage. Often, f has a single lineage
                    roll = rand(parentdist)
                    push!(edge2forest[parentedgelist[roll].number], pop!(f))
                    # other lineages: change parent probabilities in place
                    priorconc = postconc = alpha
                    parentprob = parentdist.p # parentprob !== parentdist.p, despite documentation saying so
                    for genetree in f
                        postconc += 1.0
                        probprevious = 1.0 / postconc
                        parentprob .*= priorconc * probprevious
                        parentprob[roll] += probprevious
                        roll = rand(parentdist)
                        push!(edge2forest[parentedgelist[roll].number], genetree)
                        priorconc = postconc
                    end
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

"""
    simulatecoalescent(net, nloci, nindividuals, populationsize;
        nodemapping=false, round_generationnumber=true, inheritancecorrelation=0.0)

Simulate `nloci` gene trees with `nindividuals` from each species
under the multispecies network coalescent, along network `net`,
whose branch lengths are assumed to be in **number of generations**.
`populationsize` should be a single number, assumed to be the
(haploid) effective population size Nₑ, constant across the species phylogeny.
Alternatively, `populationsize` can be a dictionary mapping the number of
each edge in `net` to its Nₑ, including an extra edge number for the
population above the root of the network.

Coalescent units are then calculated as `u=g/Nₑ` where `g` is the edge length
in `net` (number of generations), and the coalescent model is applied
using the infinite-population-size approximation.

Output: vector of gene trees with edge lengths in number of generations,
calculated as `g=uNₑ` and then rounded to be an integer, unless
`round_generationnumber` is false.

!!! warning
    When `populationsize` Nₑ is not provided as input, all edge lengths are in
    coalescent units. When `populationsize` is given as an argument, all edge
    lengths are in number of generations.
    The second method (using # generation and Nₑ as input) is a wrapper around
    the first (using coalescent units).

```jldoctest
julia> using PhyloNetworks

julia> net = readTopology("(A:500,B:500);"); # branch lengths of 100 generations

julia> Ne = Dict(e.number => 1_000 for e in net.edge);

julia> rootedgenumber = PhyloCoalSimulations.get_rootedgenumber(net)
3

julia> push!(Ne, rootedgenumber => 2_000) # Ne for population above the root
Dict{Int64, Int64} with 3 entries:
  2 => 1000
  3 => 2000
  1 => 1000

julia> using Random; Random.seed!(54321); # for replicability of example below

julia> genetrees = simulatecoalescent(net, 2, 1, Ne);

julia> writeMultiTopology(genetrees, stdout) # branch lengths: number of generations
(B:546.0,A:546.0);
(B:3155.0,A:3155.0);

```
"""
function simulatecoalescent(net::PN.HybridNetwork, nloci::Integer, nindividuals, Neff;
        nodemapping=false, round_generationnumber=true, inheritancecorrelation=0.0)
    popedgenumbers = [e.number for e in net.edge]
    push!(popedgenumbers, get_rootedgenumber(net))
    if isa(Neff, AbstractDict)
        issubset(popedgenumbers, keys(Neff)) ||
            error("populationsize is missing some edge numbers, or the root edge number.")
        valtype(Neff) <: Number || error("population sizes should be numerical")
    elseif isa(Neff, Number)
        Neff = Dict(num => Neff for num in popedgenumbers)
    else
        error("populationsize should be a number or dictionary")
    end
    # network with lengths in coalescent units
    net_coal = deepcopy(net)
    for e in net_coal.edge
        e.length = e.length / Neff[e.number]
    end
    # simulate gene trees
    genetreelist = simulatecoalescent(net_coal,nloci,nindividuals; nodemapping=true,
        inheritancecorrelation=inheritancecorrelation);
    # convert lengths to #generations in gene trees
    for gt in genetreelist
        for e in gt.edge
            len = e.length * Neff[e.inCycle]
            e.length = (round_generationnumber ? round(len) : len )
        end
        nodemapping || PN.removedegree2nodes!(gt, true) # 'true' option requires PN v0.15.1
    end
    return genetreelist
end

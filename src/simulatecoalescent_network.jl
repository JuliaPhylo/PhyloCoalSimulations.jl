"""
    simulatecoalescent(net, nloci, nindividuals=1)

Using species network `net`, simulate `nloci` gene trees
with `nindividuals` from each species.

Output: A length-`nloci` vector of gene trees

"""
function simulatecoalescent(net::PN.HybridNetwork, nloci, nindividuals=1)
    if isa(nindividuals, AbstractDict)
        # check if value types are Integers? || error("nindividuals should be integers")
        length(nindividuals) == net.numTaxa || error("nindividuals has wrong length")
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

    genetreelist = Vector{PN.HybridNetwork}(undef,nloci)
    for ilocus in 1:nloci
        for f in values(node2forest) # clean up intermediate dictionary
            empty!(f)
        end
        for f in values(edge2forest)
            empty!(f)
        end
        nextid = 1
        for nodei in nnodes:-1:1
            nn = net.nodes_changed[nodei]
            parentedgelist = parentedges[nodei] # parent population(s)
            nparents = length(parentedgelist)
            # keep track of next available node & edge ID!
            if nn.leaf
                ee = parentedgelist[1]
                f = initializetip(nn, nindividuals[nn.name], nextid)
                for e in f
                    e.inCycle = ee.number
                end
                nextid += nindividuals[nn.name]
                edge2forest[ee.number] = f
                nextid = simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
                continue 
            end
            # gather forest from children
            childedgelist = childedges[nodei] # children populations
            node2forest[nn.number] = edge2forest[childedgelist[1].number]
            f = node2forest[nn.number]
            for i in 2:length(childedgelist)
                append!(f, edge2forest[childedgelist[i].number])
            end

            if nparents == 1 # tree node but not leaf
                ee = parentedgelist[1]
                # todo: create degree-2 nodes + incomplete edges
                # name the new degree-2 nodes with string(species node number), or species node name if it exists
                # incomplete edges: set their .inCycle to ee.number
                # instead of the line below, push to edge2forest[ee.number] each new incomplete edge
                edge2forest[ee.number] = f
                nextid = simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
            elseif nparents > 1
                # partition random across all parent edges
                # todo: create degree-2 nodes + incomplete edges, name new nodes with hybrid node name
                # create the forest of new incomplete edges, to be used in place of f below
                gamma = Float64[]
                for parent in parentedgelist
                    append!(gamma, parent.gamma)
                end
                d = Distributions.Categorical(gamma)
                # below: replace f by the forest of new incomplete edges
                for genetree in f
                    roll=rand(d)
                    whereto = parentedgelist[roll] # edge in species network
                    # set genetree.inCycle to whereto.number
                    push!(edge2forest[whereto.number], genetree)
                end
                # run coalescent along each parent edge
                for ee in parentedgelist
                    nextid = simulatecoal_onepopulation!(edge2forest[ee.number], ee.length, nextid, ee.number)
                end
            else # nparents = 0: infinite root population
                if length(f) > 1 # can occur if a displayed tree's MRCA is strictly below the root, or if the network root has a single child
                    nextid = simulatecoal_onepopulation!(f, Inf, nextid)
                end
                rootnode = f[1].node[1]
                genetreelist[ilocus] = convert2tree!(rootnode)
                # may be write tree to output file?
            end
        end
    end
    return genetreelist
end

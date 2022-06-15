function simulatecoalescent(net::PN.HybridNetwork, nloci, nindividuals=1)
    if ias(nindividuals, AbstractDict)
        # check if value types are Integers? || error("nindividuals should be integers")
        length(nindividuals) == net.numTaxa || error("nindividuals has wrong length")
    elseif isa(nindividuals, Integer)
        nindividuals = Dict(n.name => nindividuals for n in net.leaf)
    else
        error("nindividuals should be an integer or vector of integers")
    end
    PN.preorder!(net)
    node2forest = Dict(n.number => PN.Edge[] for n in net.node)

    genetreelist = Vector{PN.HybridNetwork}(undef,nloci)
    nnodes = length(net.node)
    for ilocus in 1:nloci
        for f in values(node2forest) # clean up intermediate dictionary
            empty!(f)
        end
        for nodei in nnodes:-1:1
            nn = net.nodes_changed[nodei]
            parentedgelist = PN.Edge[]
            for e in nn.edge
                PN.getChild(e) === nn || continue
                push!(e, parentedgelist)
            end
            nparents = length(parentedgelist)
            nextid = 1
            # keep track of next available node & edge ID!
            if nn.leaf
                ee = parentedgelist[1]
                f = initializetip(nn, nindividuals[nn.name], nextid)
                for e in f
                    e.inCycle = ee.number
                end
                nextid += nindividuals[nn.name]
                node2forest[nn.number] = f
                simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
            elseif nparents == 1 # tree node but not leaf
                # get the index of children nodes
                # gather forests for all children populations
                # run coalescent along parent edge
            elseif nparents > 1
                # gather forest from children populations
                # partition random across all parent edges
                # run coalescent along each parent edge
            else # nparents = 0: infinite root population
                # run coalescent along infinite population
                # extract rootnode: node[1] of single edge in the forest
                genetreelist[ilocus] = convert2tree!(rootnode)
                # may be write tree to output file?
            end
        end
    end
    return genetreelist
end

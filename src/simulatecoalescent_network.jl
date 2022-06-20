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
    PN.preorder!(net)
    node2forest = Dict(n.number => PN.Edge[] for n in net.node)
    edge2forest = Dict(e.number => PN.Edge[] for e in net.edge)

    genetreelist = Vector{PN.HybridNetwork}(undef,nloci)
    nnodes = length(net.node)
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
            parentedgelist = PN.Edge[]

            #this might be a more efficient way to check edges--we only have to run the loop once
            #childedgelist = PN.Edge[]
            #for e in nn.edge
            #   PN.getChild(e) === nn ? push!(parentedgelist, e) : push!(childedgelist, e)
            #end

            for e in nn.edge
                PN.getChild(e) === nn || continue
                push!(parentedgelist, e)
            end
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
            childedgelist = PN.Edge[]
            for e in nn.edge
                PN.getChild(e) !== nn || continue
                push!(childedgelist, e)
            end

            #children = PN.getChildren(nn)
            #forest = deepcopy(node2forest[children[1].number])
            node2forest[nn.number] = edge2forest[childedgelist[1].number]
            for i in 2:length(childedgelist)
                append!(node2forest[nn.number], edge2forest[childedgelist[i].number])
            end
            
            f = node2forest[nn.number]

            if nparents == 1 # tree node but not leaf
                ee = parentedgelist[1]
                edge2forest[ee.number] = f
                nextid = simulatecoal_onepopulation!(f, ee.length, nextid, ee.number)
            elseif nparents > 1
                # partition random across all parent edges
                # (in the gene trees: do we want to make degree-2 nodes with n.inCycle = the hybrid node
                #  and new edges with e.inCycle = the appropriate hybrid parent?  not implemented right now)
                gamma = Float64[]
                for parent in parentedgelist
                    append!(gamma, parent.gamma)
                end
                d = Distributions.Categorical(gamma)
                for genetree in f
                    roll=rand(d)
                    whereto = parentedgelist[roll]
                    #parentforest = edge2forest[whereto.number]
                    #isempty(edge2forest[whereto.number]) ? edge2forest[whereto.number] = [genetree] : append!(edge2forest[whereto.number], genetree)
                    append!(edge2forest[whereto.number], [genetree])
                end
                # run coalescent along each parent edge
                for ee in parentedgelist
                    nextid = simulatecoal_onepopulation!(edge2forest[ee.number], ee.length, nextid, ee.number)
                end
            else # nparents = 0: infinite root population
                if length(f) > 1
                    nextid = simulatecoal_onepopulation!(f, Inf, nextid) # length(f)>1 needed to avoid error "there should be 2 or more lineages at the start of the infinite root population" in simulatecoal_onepopulation!
                end
                rootnode = f[1].node[1]
                genetreelist[ilocus] = convert2tree!(rootnode)
                # may be write tree to output file?
            end
        end
    end
    return genetreelist
end

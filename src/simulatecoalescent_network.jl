function simulatecoalescent(net::PN.HybridNetwork, nloci, nindividuals=1)
    if ias(nindividuals, Vector)
        eltype(nindividuals) <: Integer || error("nindividuals should be integers")
        length(nindividuals) == net.numTaxa || error("nindividuals has wrong length")
    elseif isa(nindividuals, Integer)
        nindividuals = [nindividuals for _ in 1:net.numTaxa]
    else
        error("nindividuals should be an integer or vector of integers")
    end
    PN.preorder!(net)
    # prepare data structures to store intermediate information
    # initialize tips. could those be re-used across loci?

    # write to file and/or return nloci objects?
    genetreelist = Vector{PN.HybridNetwork}(undef,nloci)
    for ilocus in 1:nloci
        rootnode = nothing # do something instead
        genetreelist[ilocus] = convert2tree!(rootnode)
    end
end

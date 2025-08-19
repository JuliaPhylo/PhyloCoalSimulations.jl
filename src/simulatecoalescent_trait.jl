"""
simulate_locuseffect_at!(node::PN.Node, transition_distribution_family, P0_distribution)
simulate_locuseffect_at!(childnode::PN.Node, transition_distribution_family, len, parenttrait)
simulate_locuseffect_below!(node::PN.Node, transition_distribution_family)
simulate_locuseffect!(Y, gtree::HybridNetwork, transition_distribution_family)

warnings:
- the nodes' field `booln5` is used to track which nodes have been simulated
  already, as a sanity check. Not used for any other purpose, only temporarily.
- no check that `Y` has the correct keys (individual names)

todo:
- argument to allow for correlated inheritance

# examples

The first example uses a Gaussian prior at the root population,
and a Brownian motion for each locus (along its own gene tree).
```jldoctest polytrait
julia> net = readnewick("((B:.5,A:.5):1.5,O:2)r;");

julia> using Distributions; root_prior = NormalCanon(0, 100) # μ=0, σ²=1/100: σ=0.1

julia> transition_distribution(x,t) = Normal(x, 0.1*sqrt(t))

julia> simulate_polygenictrait(net, 2, 1, root_prior, transition_distribution)
```

The next example uses a binary 0/1 trait at each locus, and a Markov transition
for each locus (again, along its own gene tree).
```jldoctest polytrait
julia> root_prior = Binomial(1, 0.3) # has minor allele? minor frequency 0.3

julia> function transition_distribution(x,t)
           prob_nochange = exp(-0.1*t) # substitution rate: 0.1 / coalescent unit
           prob_switch = 1 - prob_nochange
           minorprob = prob_nochange * x + (1-prob_nochange) * 0.3
           return Binomial(1, minorprob)
       end;

julia> simulate_polygenictrait(net, 1, 5, root_prior, transition_distribution)
```
"""
function simulate_polygenictrait(net::PN.HybridNetwork, args...; kwargs...)
    simulate_polygenictrait(default_rng(), net, args...; kwargs...)
end
function simulate_polygenictrait(
    rng::AbstractRNG,
    net::PN.HybridNetwork,
    nrep::Integer,
    nloci::Integer,
    P0_distribution::Distribution,
    transition_kernel;
    nindividuals=1
)
    PN.nameinternalnodes!(net, "i")
    P0_name = PN.getroot(net).name
    count(n -> n.name == P0_name, net.node) == 1 ||
        error("The root shares its name $(P0_name) with another node")
    # would cause simulate_locuseffect! to be incorrect
    net_tiplabels = PN.tiplabels(net)
    if isa(nindividuals, AbstractDict)
        issubset(net_tiplabels, keys(nindividuals)) ||
            error("nindividuals is missing some species")
        valtype(nindividuals) <: Integer ||
            error("nindividuals should be integers")
    elseif isa(nindividuals, Integer)
        nindividuals = Dict(lab => nindividuals for lab in net_tiplabels)
    else
        error("nindividuals should be an integer or dictionary of integers")
    end
    genetree_tiplabels = String[]
    for lab in net_tiplabels
        nind = nindividuals[lab]
        if nind == 1
            push!(genetree_tiplabels, lab)
        elseif nind > 1
            for i in 1:nind
                push!(genetree_tiplabels, lab * "_" * string(i))
            end
        end
    end
    lab2index = Dict(lab => i for (i,lab) in enumerate(genetree_tiplabels))
    N = length(genetree_tiplabels) # total number of individuals
    T = eltype(P0_distribution)
    p = length(P0_distribution) # number of variables
    Y = Matrix{T}(undef, N, p)  # for temp storage
    X = [zeros(typeof(zero(T)/1), N, p) for _ in 1:nrep] # X will be rescaled
    # simulate all gene trees to avoid the overhead of doing it for each rep
    locustrees = simulatecoalescent(rng, net, nloci * nrep, nindividuals;
        nodemapping=true)
    # each 'locustree' has its nodes already listed in post-order: used below
    ilocus = 0
    for irep in 1:nrep
        for _ in 1:nloci
            ilocus += 1
            tree = locustrees[ilocus]
            simulate_locuseffect!(rng, Y, lab2index, tree, P0_name,
                P0_distribution, transition_kernel)
            X[irep] .+= Y
        end
    end
    X ./= sqrt(nloci)
    return X
end

function simulate_locuseffect!(
    rng::AbstractRNG,
    Y,
    lab2index::Dict,
    gtree::PN.HybridNetwork,
    P0_name,
    P0_distribution::Distribution,
    transition_kernel,
)
    for n in gtree.node # in post-order already
        # nodes at the origin population = nodes named P0_name
        if n.name == P0_name
            simulate_locuseffect_at!(rng, Y, lab2index, n,
                P0_distribution, transition_kernel)
        end
        if n.leaf
            n.booln5 || error("the trait for $(n.name) was not simulated")
        end
    end
    return Y
end

function simulate_locuseffect_at!(
    rng::AbstractRNG,
    Y,
    lab2index::Dict,
    node::PN.Node,
    P0_distribution::Distribution,
    transition_kernel,
)
    node.booln5 &&
        error("locus effect at node number $(node.number) was already simulated")
    node.booln5 = true
    nodeY = Distributions.rand(rng, P0_distribution)
    if node.leaf
        Y[lab2index[node.name]] = nodeY
    else
        simulate_locuseffect_below!(rng, Y, lab2index, node,
            transition_kernel, nodeY)
    end
    return nothing
end
function simulate_locuseffect_at!(
    rng::AbstractRNG,
    Y,
    lab2index::Dict,
    child::PN.Node,
    transition_kernel,
    len,
    parentY,
)
    child.booln5 &&
        error("locus effect at node number $(child.number) was already simulated")
    childY = rand(rng, transition_kernel(parentY, len))
    child.booln5 = true
    if child.leaf
        Y[lab2index[child.name]] = childY
    else
        simulate_locuseffect_below!(rng, Y, lab2index, child,
            transition_kernel, childY)
    end
    return nothing
end
function simulate_locuseffect_below!(
    rng::AbstractRNG,
    Y,
    lab2index::Dict,
    node::PN.Node,
    transition_kernel,
    nodeY,
)
    for ce in node.edge
        node ≡ PN.getparent(ce) || continue # loop through child edges only
        simulate_locuseffect_at!(rng, Y, lab2index, PN.getchild(ce),
            transition_kernel, ce.length, nodeY)
    end
    return nothing
end

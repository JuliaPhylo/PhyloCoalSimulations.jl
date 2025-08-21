"""
    function simulate_polygenictrait([rng::AbstractRNG,] net, nrep, nloci,
        P0_distribution::Distribution, transition_distribution_family;
        nindividuals=1, inheritancecorrelation=0.0)
    )

Simulate a polygenetic trait controlled by multiple additive loci, each with
effect `Yi`: `X = (Y_1+...+ Y_nloci)/sqrt(nloci)`.
`nindividuals` is the number of individuals sampled from each taxon in `net`.
It can be a single integer, or a dictionary as in [`simulatecoalescent`](@ref).

Each locus effect `Y` evolves along its own gene tree, sampled from the
coalescent process, with possible correlated inheritance if
`inheritancecorrelation > 0`: see [`simulatecoalescent`](@ref).
Note that ancestors in some fixed ancestral population are
different individuals across different loci.
Branch lengths in the network are interpreted to be in coalescent units.

For each locus, `Y` at individuals in the origin population `P0`, taken to be
the root of the network, is sampled from `P0_distribution`.
Given state `y0` at the start of an edge in the gene tree,
the trait value `yt` at the end of the edge is sampled from
`transition_distribution_family(y0,t)` where `t` is the edge length,
in coalescent units.

Output:`(X, tiplabels, genetrees)` where
- `X` is a vector of matrices, of length `nrep`: one matrix for each independent
  replicate. For replicate `i`, `X[i]` is of size `N×p` where `N` is the total
  number of individuals and `p` the dimension of the trait.
- `tiplabels` lists the tips labels in gene trees, in the order in which they
  are considered in each trait matrix. That is: `X[i][j,:]` contains the traits
  in replicate `i` for individual `tiplabels[j]`.
- `genetrees` is a vector of gene trees, of length `nrep * nloci`.
  The first `nloci` trees were used for the first replicate, the next `nloci`
  trees were used for the second replicate, etc.
  If the trait is a simple scalar, then the locus effect simulated at node `n`
  in the locus tree is stored in `n.fvalue`. Other functions may erase this.

# examples

The first example uses a Gaussian prior at the root population,
and a Brownian motion for each locus (along its own gene tree).
```jldoctest polytrait
julia> net = readnewick("((B:.5,A:.5):1.5,O:2)r;");

julia> using Distributions; root_prior = NormalCanon(0, 100); # μ=0, σ²=1/100: σ=0.1

julia> transition_distribution(x,t) = Normal(x, 0.1*sqrt(t));

julia> using StableRNGs; rng = StableRNG(791); # reproducible across julia versions

julia> x,lab,gt = simulate_polygenictrait(rng, # 'rng' can be omitted
           net, 2, 1, # 2 reps, 1 locus only
           root_prior, transition_distribution);

julia> using DataFrames

julia> DataFrame(species=lab, rep1=x[1][:], rep2=x[2][:])
3×3 DataFrame
 Row │ species  rep1       rep2       
     │ String   Float64    Float64    
─────┼────────────────────────────────
   1 │ B        -0.142115  -0.0144998
   2 │ A         0.162518  -0.0736928
   3 │ O        -0.122352  -0.222733
```

The next example uses a binary 0/1 trait at each locus, and a Markov transition
for each locus (again, along its own gene tree).
```jldoctest polytrait
julia> root_prior = Binomial(1, 0.3); # has minor allele? minor frequency 0.3

julia> function transition_distribution(x,t)
           prob_nochange = exp(-0.1*t) # substitution rate: 0.1 / coalescent unit
           prob_switch = 1 - prob_nochange
           minorprob = prob_nochange * x + (1-prob_nochange) * 0.3
           return Binomial(1, minorprob)
       end;

julia> x,lab,gt = simulate_polygenictrait(rng, # rng for reproducibility only
          net, 1, 5, # 1 replicate, 5 loci: trait = (# minor alleles) / sqrt(5)
          root_prior, transition_distribution);

julia> hcat(lab, x[1])
3×2 Matrix{Any}:
 "B"  0.894427
 "A"  1.34164
 "O"  0.894427

julia> writenewick(gt[1], round=true) # gene tree for locus 1
"((O:2.0)r:1.059,(((A:0.5)i1:0.178,(B:0.5)i1:0.178):1.322)r:1.059);"
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
    nindividuals=1,
    inheritancecorrelation=0.0,
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
    storeY = size(P0_distribution) == () # 'yes' if trait is a simple scalar
    Y = Matrix{T}(undef, N, p)  # for temp storage
    X = [zeros(typeof(zero(T)/1), N, p) for _ in 1:nrep] # X will be rescaled
    # simulate all gene trees to avoid the overhead of doing it for each rep
    locustrees = simulatecoalescent(rng, net, nloci * nrep, nindividuals;
        nodemapping=true, inheritancecorrelation=inheritancecorrelation)
    # each 'locustree' has its nodes already listed in post-order: used below
    ilocus = 0
    for irep in 1:nrep
        for _ in 1:nloci
            ilocus += 1
            tree = locustrees[ilocus]
            simulate_locuseffect!(rng, Y, lab2index, tree, P0_name,
                P0_distribution, transition_kernel, storeY)
            X[irep] .+= Y
        end
    end
    X ./= sqrt(nloci)
    return X, genetree_tiplabels, locustrees
end

"""
    simulate_locuseffect!(
        rng::AbstractRNG,
        Y::Array,
        lab2index::Dict,
        gtree::PN.HybridNetwork,
        P0_name::String,
        P0_distribution::Distribution,
        transition_distribution_family
    )

Simulate one trait along its gene tree `gtree`, starting at individuals in
population named `P0_name` according to `P0_distribution`, and whose evolution,
when starting at state `y0` after time `t`, follows distribution
`transition_distribution_family(y0,t)`.

`Y` stores the simulated trait at the leaves: `Y[i]` is the trait for leaf `i`,
where `i` is some fixed index for the `leaf` node: `lab2index[leaf.name]`.

notes:
- the nodes' field `booln5` is used to track which nodes have been simulated
  already. It is only used for a sanity check, only temporarily.
- if the trait is a scalar `Float64`, then is it stored in the gene tree itself,
  at each node in its field `.fvalue`.
"""
function simulate_locuseffect!(
    rng::AbstractRNG,
    Y::Array,
    lab2index::Dict,
    gtree::PN.HybridNetwork,
    P0_name::String,
    P0_distribution::Distribution,
    transition_kernel,
    storeY::Bool,
)
    for n in gtree.node # in post-order already
        # nodes at the origin population = nodes named P0_name
        if n.name == P0_name
            simulate_locuseffect_at!(rng, Y, lab2index, n,
                P0_distribution, transition_kernel, storeY)
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
    storeY::Bool,
)
    node.booln5 &&
        error("locus effect at node number $(node.number) was already simulated")
    nodeY = Distributions.rand(rng, P0_distribution)
    if storeY
        node.fvalue = nodeY
    end
    node.booln5 = true
    if node.leaf
        Y[lab2index[node.name]] = nodeY
    else
        simulate_locuseffect_below!(rng, Y, lab2index, node,
            transition_kernel, nodeY, storeY)
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
    storeY::Bool,
)
    child.booln5 &&
        error("locus effect at node number $(child.number) was already simulated")
    childY = rand(rng, transition_kernel(parentY, len))
    if storeY
        child.fvalue = childY
    end
    child.booln5 = true
    if child.leaf
        Y[lab2index[child.name]] = childY
    else
        simulate_locuseffect_below!(rng, Y, lab2index, child,
            transition_kernel, childY, storeY)
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
    storeY::Bool,
)
    for ce in node.edge
        node ≡ PN.getparent(ce) || continue # loop through child edges only
        simulate_locuseffect_at!(rng, Y, lab2index, PN.getchild(ce),
            transition_kernel, ce.length, nodeY, storeY)
    end
    return nothing
end

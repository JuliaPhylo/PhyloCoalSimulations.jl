module PhyloCoalSimulations

using Distributions
import PhyloNetworks     # brings PhyloNetworks in scope, but none of its names
import Random: AbstractRNG, default_rng

const PN = PhyloNetworks # import PhyloNetworks as PN requires julia v1.6

export
simulatecoalescent,
population_mappedto,
gene_edgemapping!

include("simulatecoalescent_onepop.jl")
include("simulatecoalescent_network.jl")
include("utils.jl")

end

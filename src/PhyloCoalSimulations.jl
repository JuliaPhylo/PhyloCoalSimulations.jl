module PhyloCoalSimulations

using Distributions
import PhyloNetworks     # brings PhyloNetworks in scope, but none of its names
const PN = PhyloNetworks # import PhyloNetworks as PN requires julia v1.6

export
simulatecoalescent

include("simulatecoalescent_onepop.jl")
include("simulatecoalescent_network.jl")

end

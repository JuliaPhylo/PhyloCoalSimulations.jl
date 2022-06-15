module PhyloCoalSimulations

using Distributions
import PhyloNetworks as PN # brings PN in scope, but none of its names

export
simulatecoalescent

include("simulatecoalescent_onepop.jl")
include("simulatecoalescent_network.jl")

end

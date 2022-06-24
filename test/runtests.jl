using PhyloCoalSimulations
import PhyloCoalSimulations as PCS
using Test
using Random
using PhyloNetworks
import PhyloNetworks as PN

@testset "PhyloCoalSimulations.jl" begin
    include("test_onepopulation.jl")
    include("test_multispeciesnetwork.jl")
end

using PhyloCoalSimulations
using Test

@testset "PhyloCoalSimulations.jl" begin
    # tests could be broken across files
    @test simulatecoal_onepopulation(1, 2.0) == 1
    @test isnothing(simulatecoal_onepopulation(0, 2.0))
end

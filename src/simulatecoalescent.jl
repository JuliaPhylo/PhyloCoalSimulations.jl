"""
    simulatecoal_onepopulation(nlineage, population_length)

write description here

# examples

```jldoctest
julia> simulatecoal_onepopulation(1, 2.0)
hello population: 1 to start with
1

julia> simulatecoal_onepopulation(0, 2.0)

```
"""
function simulatecoal_onepopulation(nlineage::Int, poplen::AbstractFloat)
    nlineage < 1 && return(nothing)
    println("hello population: $nlineage to start with")
    return nlineage
end

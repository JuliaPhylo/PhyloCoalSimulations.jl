"""
    population_mappedto(edge or node)

Identifier of the population (edge in the species network) that a gene tree's
edge or a node is mapped onto, or `nothing` if not mapped. For example,
coalescent nodes in gene trees map to a *node* in the species phylogeny, instead
of mapping to an *edge*.
"""
population_mappedto(e::Union{PN.Edge,PN.Node}) = (e.inte1 == -1 ? nothing : e.inte1)

"""
    ismappingnode(node)

Boolean: true if `node` is of degree 2, has a single child, and has a name.
(The root is of degree-2 but is not a mapping node).
"""
ismappingnode(n::PN.Node) = length(n.edge) == 2 && PN.hassinglechild(n) && n.name != ""

"""
    mappingnodes(gene tree)

Type to define an iterator over degree-2 mapping nodes in a gene tree, assuming these
degree-2 nodes (other than the root) have a name to map them to nodes in a
species phylogeny. See [`ismappingnode`](@ref PhyloCoalSimulations.mappingnodes).
"""
struct mappingnodes
    phy::PN.HybridNetwork
end
function Base.iterate(mn::mappingnodes, state=1)
    next = iterate(mn.phy.node, state)
    while next !== nothing
        (n, st) = next
        ismappingnode(n) && break
        next = iterate(mn.phy.node, st)
    end
    return next
end
Base.IteratorSize(::Type{mappingnodes}) = Base.SizeUnknown()
Base.eltype(::Type{mappingnodes}) = PN.Node
# see https://docs.julialang.org/en/v1/manual/interfaces/ for interators

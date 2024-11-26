"""
    population_mappedto(edge or node)

Identifier of the population (edge in the species network) that a gene tree's
edge or a node is mapped onto, or `nothing` if not mapped. For example,
coalescent nodes in gene trees map to a *node* in the species phylogeny, instead
of mapping to an *edge*.
"""
population_mappedto(e::PN.Node) = (e.intn1 == -1 ? nothing : e.intn1)
population_mappedto(e::PN.Edge) = (e.inte1 == -1 ? nothing : e.inte1)

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


"""
    gene_edgemapping!(gene_tree, species_network, checknames=true)

Given a gene tree with labeled internal nodes that map to a species phylogeny
(a species tree or a species network),
this function maps each gene edge to the species edge that it is contained "within".
Gene edge mappings are stored in the `inte1` field of each gene tree edge,
but it's best to access this mapping via [`population_mappedto`](@ref).

Assumption: the `species_network` has unique node names to uniquely identify
the speciation and reticulation events; and the `gene_tree` has degree-2 nodes
with matching names, to indicate which species event each degree-2 node
corresponds to.

The `checknames` argument takes a boolean, and, if `true`, 
then the function will check that both the species and the gene phylogeny
have the same internal node names. The mapping of edges is recovered from
matching names between nodes in the gene tree and nodes in the species network.
"""
function gene_edgemapping!(
    gene::PN.HybridNetwork,
    species::PN.HybridNetwork,
    checknames::Bool=true
)
    # get all named non-leaf nodes in the gene tree: nodes to traverse for edge groups
    search_nds = gene.node[(x->x.name != "" && !x.leaf).(gene.node)]
    if checknames
        gene_nodenames = (x->x.name).(gene.node)
        species_nodenames = push!((x->x.name).(species.node),"")
        species_nodename_set = Set(species_nodenames)
        length(species_nodenames) == length(species_nodename_set) ||
            error("The species phylogeny does not a unique name for each node")
        issetequal(species_nodename_set, gene_nodenames) ||
            error("The gene and species phylogeny have different sets of node names")
    end
    # if the root has no name, also add the root as a starting point
    generoot = PN.getroot(gene)
    generoot.name == "" && push!(search_nds, generoot)

    # dictionary that tracks all edges found within a given species edge
    # Species edges are defined by the named parent and child nodes
    edge_groups = Dict{Tuple{String,String},Vector{PhyloNetworks.Edge}}()
    edge_group_sp_names= (x->(PN.getparent(x).name,PN.getchild(x).name)).(species.edge)
    # add an edge group for the 'root edge' to handle coalescences above the species root
    push!(edge_group_sp_names,(generoot.name, PN.getroot(species).name))
    # initialize edge groups with empty vectors
    (x -> edge_groups[x]=Vector{PhyloNetworks.Edge}()).(edge_group_sp_names)

    for search_nd in search_nds 
        nd_group = PN.getchildren(search_nd) # The nodes found within an edge group
        edge_group = PN.Edge[]
        chld_name = nothing # The name of the child node of the species edge
        for nd in nd_group
            push!(edge_group,PN.getparentedge(nd))
            if nd.name==""
               append!(nd_group,PN.getchildren(nd)) 
            else
                if isnothing(chld_name)
                    chld_name=nd.name
                else
                   chld_name != nd.name && error("Found two different child names of the species edge within an edge group.") 
                end
            end
        end
        append!(edge_groups[(search_nd.name,chld_name)],edge_group)
    end

    inte1_nos = (x-> x.number ).(species.edge)
    inte1_nos = push!(inte1_nos,length(species.edge)+1)
    for (number,e_names) in zip(inte1_nos,edge_group_sp_names)
        gene_edges = edge_groups[e_names]
        (x-> x.inte1 = number).(gene_edges)
    end
    return nothing
end

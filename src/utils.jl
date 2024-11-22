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


"""
encode_edges!(gene tree, species tree)

Given a gene tree with labeled internal nodes that map to a species phylogeny,
this function maps each gene edge to the species edge that it is contained "within".
Gene edge mappings are stored in the `inCycle` field of each gene tree `Edge`.

The `checknames` argument takes a boolean, and, if `true`, 
then the function will check that both the species and the gene phylogeny
have the same internal node names to create a mapping.
"""

function encode_edges!(gene::PN.HybridNetwork,species::PN.HybridNetwork,checknames=true::Bool)
    #Get all named non-leaf nodes on the gene tree.
    search_nds = gene.node[(x->x.name != "" && !x.leaf).(gene.node)] ##These are the  nodes that we will traverse for edge groups
    if checknames
        gene_nodenames = (x->x.name).(gene.node)
        species_nodenames = push!((x->x.name).(species.node),"")
        !issetequal(species_nodenames,gene_nodenames) && error("The gene and species phylogeny have different sets of node names")
    end ##else assume only node corresponding to the species tree are named in the gene tree
    gene.node[gene.root].name == "" && push!(search_nds,gene.node[gene.root]) ## if not named, also add the root as a starting point

    # dictionary that tracks all edges found within a given species edge
    # Species edges are defined by the named parent and child nodes
    edge_groups = Dict{Tuple{String,String},Vector{PhyloNetworks.Edge}}()
    edge_group_sp_names= (x->(PN.getparent(x).name,PN.getchild(x).name)).(species.edge)
    push!(edge_group_sp_names,(gene.node[gene.root].name,species.node[species.root].name)) ## Add an edge group for the 'root edge' to handle coalescences beyond the species root. 
    (x -> edge_groups[x]=Vector{PhyloNetworks.Edge}()).(edge_group_sp_names) # initialize edge groups with empty vectors. 

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

    inCycle_nos = (x-> x.number ).(species.edge)
    inCycle_nos = push!(inCycle_nos,length(species.edge)+1)
    for (number,e_names) in zip(inCycle_nos,edge_group_sp_names)
        gene_edges = edge_groups[e_names]
        (x-> x.inCycle = number).(gene_edges)
    end
    return nothing
end
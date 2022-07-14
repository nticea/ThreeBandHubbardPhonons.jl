using Graphs
using GraphRecipes
using Plots 
include("cuprate_lattice.jl")

"""
Build graph from lattice. This can be useful for analyzing the bond structure
with graph algorithms and visualization of the lattice as a sanity check
"""
function build_graph(lattice::LatticeCuprate)::SimpleGraph
    edges = map(bond -> Edge(bond.s1, bond.s2), lattice)
    return SimpleGraph(edges)
end

Coord = Tuple{<:Real, <:Real}

"""
Builds a dictionary mapping the site index to cartesian coordinates
"""
function get_coords(lattice::LatticeCuprate)
    isempty(lattice) && return []

    # make index => (x, y) pairs
    function getcoord(bond)
        if bond.type=="dpx"
            x1 = 2*(bond.x1-1)
            y1 = 2*(bond.y1-1)
            x2 = 2*(bond.x2-1) + 1
            y2 = 2*(bond.y2-1)
            return [bond.s1 => (x1, y1), bond.s2 => (x2, y2)]
        elseif bond.type=="dpy"
            x1 = 2*(bond.x1-1)
            y1 = 2*(bond.y1-1)
            x2 = 2*(bond.x2-1) 
            y2 = 2*(bond.y2-1) + 1 
            return [bond.s1 => (x1, y1), bond.s2 => (x2, y2)]
        elseif bond.type=="pxpy"
            x1 = 2*(bond.x1-1) + 1
            y1 = 2*(bond.y1-1) 
            x2 = 2*(bond.x2-1)
            y2 = 2*(bond.y2-1) + 1
            return [bond.s1 => (x1, y1), bond.s2 => (x2, y2)]
        else
            @error "Bond type is not recognized"
        end
        
    end
    allcoords = vcat(getcoord.(lattice)...)

    # create dictionary to remove duplicates
    dict = Dict(allcoords)

    # create array or coordinates indexed by the site index
    indices = keys(dict)
    maxidx = maximum(indices)
    coordlist = Vector{Coord}(undef, maxidx)
    for i in 1:maxidx
        coordlist[i] = get(dict, i, (-1, -1))
    end

    x = [coord[1] for coord in coordlist]
    y = [coord[2] for coord in coordlist]

    return x,y
end

function get_defaults(lattice::LatticeCuprate)
    if lattice[1].type == "dpx" || lattice[1].type == "dpy"
        return (
            curves=false,
            curvature_scalar=0.2,
            nodesize=1,
            nodeshape=:circle,
            nodecolor=:purple,
        )
    elseif lattice[1].type == "pxpy" 
        return (
            curves=true,
            curvature_scalar=-0.2,
            nodesize=1,
            nodeshape=:rect,
            nodecolor=:orange,
        )
    end
end

function get_sign(lattice::LatticeCuprate)
    getsign(bond) = bond.sign
    return getsign.(lattice)
end

"""
Visualize ITensor lattice
Keyword arguments are GraphRecipes arguments: https://docs.juliaplots.org/stable/generated/graph_attributes/
"""
function visualize(lattice::LatticeCuprate; kwargs...)
    if isempty(lattice)
        @warn "Empty lattice — nothing to visualize"
        return
    end

    sign = get_sign(lattice)
    lattice_pos = lattice[sign .> 0]
    lattice_neg = lattice[sign .< 0]

    defaultargs = get_defaults(lattice_pos)
    graph = build_graph(lattice_pos)
    x,y = get_coords(lattice_pos)
    graphplot(graph; x, y, edgecolor=:red, defaultargs..., kwargs...)

    defaultargs = get_defaults(lattice_neg)
    graph = build_graph(lattice_neg)
    x,y = get_coords(lattice_neg)
    graphplot!(graph; x, y, edgecolor=:blue, defaultargs..., kwargs...)
end

function visualize!(lattice::LatticeCuprate; kwargs...)
    if isempty(lattice)
        @warn "Empty lattice — nothing to visualize"
        return
    end

    sign = get_sign(lattice)
    lattice_pos = lattice[sign .> 0]
    lattice_neg = lattice[sign .< 0]

    defaultargs = get_defaults(lattice_pos)
    graph = build_graph(lattice_pos)
    x,y = get_coords(lattice_pos)
    graphplot!(graph; x, y, edgecolor=:red, defaultargs..., kwargs...)

    defaultargs = get_defaults(lattice_neg)
    graph = build_graph(lattice_neg)
    x,y = get_coords(lattice_neg)
    graphplot!(graph; x, y, edgecolor=:blue, defaultargs..., kwargs...)
end

# function visualize!(lattice::LatticeCuprate; use_lattice_coords=true, kwargs...)
#     if isempty(lattice)
#         @warn "Empty lattice — nothing to visualize"
#         return
#     end

#     defaultargs = get_defaults(lattice)

#     graph = build_graph(lattice)

#     if use_lattice_coords
#         coords = get_coords(lattice)
#         x = [coord[1] for coord in coords]
#         y = [coord[2] for coord in coords]
#         graphplot!(graph; x, y, defaultargs..., kwargs...)
#     else
#         graphplot!(graph; defaultargs..., kwargs...)
#     end
# end
"""
    abstract type AbstractMesh{N,T}

An abstract mesh structure in dimension `N` with primite data of type `T` (e.g.
`Float64` for double precision representation).

Concrete subtypes of `AbstractMesh` should implement [`ElementIterator`](@ref)
for accessing the mesh elements.

See also: [`LagrangeMesh`](@ref)
"""
abstract type AbstractMesh{N,T} end

ambient_dimension(::AbstractMesh{N}) where {N} = N

function Base.show(io::IO, msh::AbstractMesh)
    print(io, "$(typeof(msh)) containing:")
    for E in element_types(msh)
        iter = elements(msh, E)
        print(io, "\n\t $(length(iter)) elements of type ", E)
    end
    return io
end

"""
    struct ElementIterator{E,M}

Iterator for all elements of type `E` on a mesh of type `M`.

Besides the methods listed in the [iterator
iterface](https://docs.julialang.org/en/v1/manual/interfaces/) of `Julia`, some
functions also require the `getindex(iter,i::Int)` method for accessing the
`i`-th element directly.
"""
struct ElementIterator{E,M}
    mesh::M
end

"""
    elements(msh::AbstractMesh,E::DataType)

Return an iterator for all elements of type `E` on a mesh `msh`.
"""
elements(msh::AbstractMesh, E::DataType) = ElementIterator{E,typeof(msh)}(msh)

"""
    struct LagrangeMesh{N,T} <: AbstractMesh{N,T}

Data structure representing a generic mesh in an ambient space of dimension `N`,
with data of type `T`.

The `LagrangeMesh` can, in principle, store elements of any type. Those are given
as a key in the `elements` dictionary, and the value is a data structure which
is capable of reconstructing the elements. For example, for a Lagrange element
described by `p` nodes, the value is a `p×Nel` matrix of integer, where each
columns is a list of tags for the nodes of the element. The nodes are stored in
the `nodes` field.
"""
struct LagrangeMesh{N,T} <: AbstractMesh{N,T}
    nodes::Vector{SVector{N,T}}
    # for each element type (key), return the connectivity matrix
    etype2mat::Dict{DataType,Matrix{Int}}
    # mapping from entity to a dict containing (etype=>tags)
    ent2etags::Dict{EntityKey,Dict{DataType,Vector{Int}}}
end

# empty constructor
function LagrangeMesh{N,T}() where {N,T}
    return LagrangeMesh{N,T}(
        SVector{N,T}[],
        Dict{DataType,Matrix{Int}}(),
        Dict{EntityKey,Dict{DataType,Vector{Int}}}(),
    )
end

"""
    element_types(msh::AbstractMesh)

Return the element types present in the `msh`.
"""
element_types(msh::LagrangeMesh) = keys(msh.etype2mat)

nodes(msh::LagrangeMesh) = msh.nodes
ent2etags(msh::LagrangeMesh) = msh.ent2etags
etype2mat(msh::LagrangeMesh) = msh.etype2mat

function ent2nodetags(msh::LagrangeMesh, ent::EntityKey)
    tags = Int[]
    for (E, t) in msh.ent2etags[ent]
        append!(tags, view(msh.etype2mat[E], :, t))
    end
    return tags
end

entities(msh::LagrangeMesh) = keys(msh.ent2etags)

"""
    domain(msh::LagrangeMesh)

Set of all entities covered by the mesh.
"""
domain(msh::LagrangeMesh) = Domain(entities(msh))

"""
    dom2elt(m::LagrangeMesh,Ω,E)::Vector{Int}

Compute the element indices `idxs` of the elements of type `E` composing `Ω`.
"""
function dom2elt(m::LagrangeMesh, Ω::Domain, E::DataType)
    idxs = Int[]
    for ent in entities(Ω)
        tags = get(m.ent2etags[ent], E, Int[])
        append!(idxs, tags)
    end
    return idxs
end

# implement the interface for ElementIterator of lagrange elements on a generic
# mesh. The elements are constructed on the flight based on the global nodes and
# the connectivity list stored
function Base.length(iter::ElementIterator{E,<:LagrangeMesh}) where {E}
    tags = iter.mesh.etype2mat[E]::Matrix{Int}
    _, Nel = size(tags)
    return Nel
end

function Base.getindex(iter::ElementIterator{E,<:LagrangeMesh}, i::Int) where {E}
    tags = iter.mesh.etype2mat[E]::Matrix{Int}
    node_tags = view(tags, :, i)
    vtx = view(iter.mesh.nodes, node_tags)
    el = E(vtx)
    return el
end

function Base.iterate(iter::ElementIterator{<:LagrangeElement,<:LagrangeMesh}, state = 1)
    state > length(iter) && (return nothing)
    return iter[state], state + 1
end

# convert a mesh to 2d by ignoring third component. Note that this also requires
# converting various element types to their 2d counterpart. These are needed
# because some meshers like gmsh always create three-dimensional objects, so we
# must convert after importing the mesh
function _convert_to_2d(mesh::LagrangeMesh{3,T}) where {T}
    # create new dictionaries for elements and ent2etagsdict with 2d elements as keys
    new_etype2mat = empty(mesh.etype2mat)
    new_ent2etags = empty(mesh.ent2etags)
    for (E, tags) in mesh.etype2mat
        E2d = _convert_to_2d(E)
        new_etype2mat[E2d] = tags
    end
    for (ent, dict) in mesh.ent2etags
        new_dict = empty(dict)
        for (E, tags) in dict
            E2d = _convert_to_2d(E)
            new_dict[E2d] = tags
        end
        new_ent2etags[ent] = new_dict
    end
    # construct new 2d mesh
    return LagrangeMesh{2,T}([x[1:2] for x in mesh.nodes], new_etype2mat, new_ent2etags)
end

function _convert_to_2d(::Type{LagrangeElement{R,N,SVector{3,T}}}) where {R,N,T}
    return LagrangeElement{R,N,SVector{2,T}}
end
_convert_to_2d(::Type{SVector{3,T}}) where {T} = SVector{2,T}

"""
    struct SubMesh{N,T} <: AbstractMesh{N,T}

View into a `parent` mesh over a given `domain`.

A submesh implements the interface for `AbstractMesh`; therefore you can iterate
over elements of the submesh just like you would with a mesh.

Construct `SubMesh`s using `view(parent,Ω::Domain)`.
"""
struct SubMesh{N,T} <: AbstractMesh{N,T}
    parent::LagrangeMesh{N,T}
    domain::Domain
    # etype2etags maps E => indices of elements in parent mesh contained in the
    # submesh
    etype2etags::Dict{DataType,Vector{Int}}
    function SubMesh(mesh::LagrangeMesh{N,T}, Ω::Domain) where {N,T}
        etype2etags = Dict{DataType,Vector{Int}}()
        for E in element_types(mesh)
            # add the indices of the elements of type E in the submesh. Skip if
            # empty
            idxs = dom2elt(mesh, Ω, E)
            isempty(idxs) || (etype2etags[E] = idxs)
        end
        return new{N,T}(mesh, Ω, etype2etags)
    end
end
Base.view(m::LagrangeMesh, Ω::Domain) = SubMesh(m, Ω)
Base.view(m::LagrangeMesh, ent::EntityKey) = SubMesh(m, Domain(ent))

Base.collect(msh::SubMesh) = msh.parent[msh.domain]

ambient_dimension(::SubMesh{N}) where {N} = N

geometric_dimension(msh::AbstractMesh) = geometric_dimension(domain(msh))

domain(msh::SubMesh) = msh.domain
entities(msh::SubMesh) = entities(domain(msh))

element_types(msh::SubMesh) = keys(msh.etype2etags)

ent2nodetags(msh::SubMesh, ent::EntityKey) = ent2nodetags(msh.parent, ent)

"""
    nodetags(msh::SubMesh)

Return the tags of the nodes in the parent mesh belonging to the submesh.
"""
function nodetags(msh::SubMesh)
    tags = Int[]
    for ent in entities(msh)
        append!(tags, ent2nodetags(msh, ent))
    end
    return unique!(tags)
end
nodetags(msh::LagrangeMesh) = collect(1:length(msh.nodes))

"""
    nodes(msh::SubMesh)

A view of the nodes of the parent mesh belonging to the submesh. The ordering is
given by the [`nodetags`](@ref) function.
"""
function nodes(msh::SubMesh)
    tags = nodetags(msh)
    return view(msh.parent.nodes, tags)
end

# ElementIterator for submesh
function Base.length(iter::ElementIterator{E,<:SubMesh}) where {E<:LagrangeElement}
    submesh = iter.mesh
    idxs    = submesh.etype2etags[E]::Vector{Int}
    return length(idxs)
end

function Base.getindex(
    iter::ElementIterator{E,<:SubMesh},
    i::Int,
) where {E<:LagrangeElement}
    submsh = iter.mesh
    p_msh  = submsh.parent # parent mesh
    idxs   = submsh.etype2etags[E]::Vector{Int}
    iglob  = idxs[i] # global index of element in parent mesh
    iter   = elements(p_msh, E) # iterator over parent mesh
    return iter[iglob]
end

function Base.iterate(iter::ElementIterator{<:LagrangeElement,<:SubMesh}, state = 1)
    state > length(iter) && (return nothing)
    return iter[state], state + 1
end

function Base.getindex(msh::LagrangeMesh, Ω::Domain)
    nodes = empty(msh.nodes)
    etype2mat = empty(msh.etype2mat)
    ent2etags = empty(msh.ent2etags)
    foreach(ent -> ent2etags[ent] = Dict{DataType,Vector{Int}}(), entities(Ω))
    glob2loc = Dict{Int,Int}()
    for E in element_types(msh)
        connect = msh.etype2mat[E]::Matrix{Int}
        np, _ = size(connect)
        mat = Int[]
        for ent in entities(Ω)
            etags = Int[]
            haskey(msh.ent2etags[ent], E) || continue
            for (iloc, i) in enumerate(msh.ent2etags[ent][E])
                push!(etags, iloc)
                for j in view(connect, :, i)
                    if !haskey(glob2loc, j) # new node
                        push!(nodes, msh.nodes[j])
                        glob2loc[j] = length(nodes)
                    end
                    push!(mat, glob2loc[j]) # push local index of node
                end
            end
            push!(ent2etags[ent], E => etags)
        end
        isempty(mat) || (etype2mat[E] = reshape(mat, np, :))
    end
    return LagrangeMesh(nodes, etype2mat, ent2etags)
end

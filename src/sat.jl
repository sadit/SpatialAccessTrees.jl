# this file is part of SpatialAccessTrees.jl

export Sat, SatInitialPartition, RandomInitialPartition,
    RandomSortSat, ProximalSortSat, DistalSortSat,
    search, searchbatch, index!, satpermutation, satpermutation!, permute
using Polyester

abstract type AbstractSortSat end
struct RandomSortSat <: AbstractSortSat end
struct ProximalSortSat <: AbstractSortSat end
struct DistalSortSat <: AbstractSortSat end

abstract type AbstractInitialPartition end

struct SatInitialPartition <: AbstractInitialPartition end

struct RandomInitialPartition <: AbstractInitialPartition
    nparts::Int
    shuffle::Bool
    RandomInitialPartition(; nparts=max(4, Threads.nthreads()), shuffle=false) = new(nparts, shuffle)
end

"""
    struct Sat{DT<:SemiMetric,DBT<:AbstractDatabase} <: AbstractSearchIndex
        dist::DT
        db::DBT
        root::UInt32
        # parents::Vector{UInt32}
        children::Vector{Union{Nothing,Vector{UInt32}}}
        cov::Vector{Float32} # leafs: ``- d(parent, leaf)``, internal: ``\\max \\{d(parent, u) | u \\in children(parent)\\}``
    end

Spatial Access Tree data structure. Please see `Sat` function for high level constructors
"""
struct Sat{DT<:SemiMetric,DBT<:AbstractDatabase} <: AbstractSearchIndex
    dist::DT
    db::DBT
    root::UInt32
    # parents::Vector{UInt32}
    children::Vector{Union{Nothing,Vector{UInt32}}}
    cov::Vector{Float32} # leafs: d(parent, leaf), internal: max {d(parent, u) | u \in children(parent)}
end

function Base.show(io::IO, sat::Sat)
    io = IOContext(io, :compact => true, :limit => true)
    println(io, "{\n  ", typeof(sat), " n=", length(sat), "\n")

    for j in eachindex(sat.children)
        print(io, "  [", j, "] ")
        show(io, sat.children[j])
        print(io, "\n")
        if j >= 30
            println(io, "...")
            break
        end
    end

    println(io, sat.cov)
    println(io, "}")
end

function Sat(
        sat::Sat;
        dist=sat.dist,
        db=sat.db,
        root=sat.root,
        # parents=sat.parents,
        children=sat.children,
        cov=sat.cov
    )

    # Sat(dist, db, convert(UInt32, root), parents, children, cov)
    Sat(dist, db, convert(UInt32, root), children, cov)
end

"""
    Sat(db::AbstractDatabase; dist::SemiMetric=L2Distance(), root=1)

Prepares the metric data structure. After calling this constructor, please call `index!`.

# Arguments

- `db`: database to index

# Keyword arguments
- `dist`: distance function, defaults to L2Distance()
- `root`: The dataset's element to be used as root
"""
function Sat(db::AbstractDatabase; dist::SemiMetric=L2Distance(), root=1)
    n = length(db)
    # P = zeros(UInt32, n)
    C = Union{Nothing,Vector{UInt32}}[nothing for _ in 1:n]
    cov = Vector{Float32}(undef, n)
    Sat(dist, db, convert(UInt32, root), C, cov)
end

@inline getpools(::Sat) = nothing
@inline database(sat::Sat) = sat.db
@inline database(sat::Sat, i) = sat.db[i]
@inline distance(sat::Sat) = sat.dist
@inline Base.length(sat::Sat) = length(sat.cov)

"""
    index!(sat::Sat; <kwargs...>)
    index!(sat::Sat, ipart::SatInitialPartition; <kwargs...>)
    index!(sat::Sat, ipart::RandomInitialPartition; <kwargs...>)

Performs the indexing of the referenced dataset in the tree. It supports limited forms of multithreading, induced by initial partitioning schemes.

# Arguments
- `sat`: The metric data structure.
- `ipart`: initial partitioning scheme for the tree. It supports the following kinds of objects:
    - `SatInitialPartition()`: Traditional construction, default value. Each part is a SAT partition and will be processed by a thread on multithreading setups.
    - `RandomInitialPartition(nparts=Threads.nthreads(), shuffle=false)`:
        construction that divides the dataset (randomly if `suffle=true`) in `nparts` disjoint parts. The resulting structure violates the SAT partitioning in a whole
        and creates a kind of SAT forest that are fine SAT partitions. Useful to limit the height of the tree and for multiprocessing purporses, i.e., each part will be processed by a thread.

# Keyword arguments
- `sortsat`: The strategy to create the spatial access tree, it heavily depends on the order of elements while it is build. It accepts:
   - `RandomSortSat()`: children are randomized (default value)
   - `ProximalSortSat()`: classical approach, near elements are put first.
   - `DistalSortSat()`: recent approach, distant elements are put first.
- `minleaf`: Minimum number of children to perform a spatial access separation (half space partitioning)
"""
function index!(
        sat::Sat, ipart::SatInitialPartition=SatInitialPartition();
        sortsat::AbstractSortSat=RandomSortSat(),
        minleaf::Int=log2(ceil(database(sat)))
    )
    n = length(sat)
    p::UInt32 = sat.root
    sat.children[p] = collect(UInt32, Iterators.flatten((1:p-1, p+1:n)))

    if Threads.nthreads() == 1
        D = Vector{Tuple{Float32,UInt32}}(undef, n)
        index_loop!(sat, sortsat, D, minleaf, UInt32[p])
    else
        let D = Vector{Tuple{Float32,UInt32}}(undef, n)
            index_sat_neighbors!(sat, sortsat, sat.children[p], p, D, 0)
        end
        C = sat.children[p]

        @batch per=thread minbatch=1 for i in eachindex(C)
            c = C[i]
            if sat.children[c] !== nothing
                D = Vector{Tuple{Float32,UInt32}}(undef, length(sat.children[c]))
                index_loop!(sat, sortsat, D, minleaf, UInt32[c])
            end
        end
    end

    sat
end

function index!(
        sat::Sat,
        ipart::RandomInitialPartition;
        sortsat::AbstractSortSat=RandomSortSat(),
        minleaf::Int=log2(ceil(database(sat)))
    )
    n = length(sat)
    nparts = 8ipart.nparts > n ? ceil(Int, ipart.nparts / 8) : ipart.nparts
    nparts == 1 && return index!(sat, SatInitialPartition(); sortsat, minleaf)

    p::UInt32 = sat.root    
    P = collect(UInt32, Iterators.flatten((1:p-1, p+1:n)))
    n = length(P)
    ipart.shuffle && shuffle!(P)
    sat.children[p] = P[1:nparts]
    m = ceil(Int, (n - nparts) / nparts)

    @batch per=thread minbatch=1 for i in 1:nparts
        sp = nparts + (i-1) * m
        c = P[i]
        C = sat.children[c] = P[sp+1:min(n, sp+m)]
        D = Vector{Tuple{Float32,UInt32}}(undef, length(C))
        index_loop!(sat, sortsat, D, minleaf, UInt32[c])
    end

    cov = sat.cov
    cov[p] = 0f0
    for c in (sat.children[p])::Vector{UInt32}
        cov[p] = max(cov[p], abs(cov[c]))
    end

    sat
end

function index_loop!(sat::Sat, sortsat::AbstractSortSat, D::AbstractVector, minleaf::Int, queue::Vector{UInt32})
    while length(queue) > 0
        p = pop!(queue)
        index_sat_neighbors!(sat, sortsat, sat.children[p], p, D, minleaf)
        # @assert sat.children[p] !== nothing
        @inbounds for c in (sat.children[p])::Vector{UInt32}
            C = sat.children[c]
            C !== nothing && push!(queue, c)
        end
    end

    sat
end

function index_sat_neighbors!(sat::Sat, sortsat::AbstractSortSat, C::Nothing, p::UInt32, D::Vector, minleaf::Integer)
    # do nothing
end

function index_sat_neighbors!(sat::Sat, sortsat::AbstractSortSat, C::AbstractVector, p::UInt32, D::Vector, minleaf::Integer)
    # note: D is a cache of distances and objects, it is used in two ways in this function
    n = length(C)

    resize!(D, n)
    parent = database(sat, p)
    dist = distance(sat)

    # computing distance to its parent (stored in D)
    sat.cov[p] = 0f0

    for i in eachindex(C) #(i, c) in enumerate()
        c = C[i]
        d = evaluate(dist, parent, database(sat, c))
        D[i] = (convert(Float32, d), c)
        # not thread-safe:
        sat.cov[p] = max(sat.cov[p], d) # covering radius
    end

    T = typeof(sortsat)
    T === RandomSortSat && minleaf == 0 || sort!(D, by=first)
    minleaf = min(n, minleaf)

    if minleaf > 0
        @inbounds for i in 1:minleaf  # mandatory leafs
            (d_, i_) = D[i]
            D[i] = (-d_, i_)
        end
    end

    if T === RandomSortSat
        shuffle!(D)
    elseif T === DistalSortSat
        reverse!(D)
    end

    # computing nearest neighbors of $child \in D$ (using previous D and storing the new set on D)
    empty!(C)

    for (d_, i_) in D
        d_ <= 0  && continue # negative distances encode mandatory leafs, see next outside-for-loop

        res = getknnresult(1)
        push_item!(res, IdWeight(p, d_))  # insert parent
        child = database(sat, i_)

        for j in C
            d = evaluate(dist, child, database(sat, j))
            push_item!(res, IdWeight(j, d))
        end

        nn = argmin(res)
        sat.cov[i_] = minimum(res) # distance to parent, marked as negative
        # sat.parents[i_] = nn
        if sat.children[nn] === nothing
            sat.children[nn] = UInt32[i_]
        else
            push!(sat.children[nn], i_)
        end
    end

    for (d_, i_) in D
        if d_ <= 0  # negative distances encode mandatory leafs
            sat.cov[i_] = abs(d_)
            push!(C, i_)
            continue 
        end
    end
end

#=
isleaf(sat::Sat, i::Integer) = sat.cov[i] < 0
isinner(sat::Sat, i::Integer) = sat.cov[i] >= 0
isroot(sat::Sat, i::Integer) = sat.root == i
=#
function searchtree(sat::Sat, q, p::Integer, res::KnnResult)
    cost = 1
    dist = distance(sat)
    dqp = evaluate(dist, q, database(sat, p))
    push_item!(res, IdWeight(p, dqp))

    if sat.children[p] !== nothing # inner node
        if length(res) < maxlength(res) || dqp < maximum(res) + sat.cov[p]
            C = sat.children[p]::Vector{UInt32}
            for c in C
                cost += searchtree(sat, q, c, res)
            end
        end
    end

    cost
end

function search(sat::Sat, q::T, res::KnnResult; pools=nothing) where T
    cost = searchtree(sat, q, sat.root, res)
    (; res, cost)
end


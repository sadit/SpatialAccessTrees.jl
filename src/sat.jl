# this file is part of SpatialAccessTrees.jl

export Sat, RandomSortSat, ProximalSortSat, DistalSortSat, search, searchbatch, index!

abstract type AbstractSortSat end
struct RandomSortSat <: AbstractSortSat end
struct ProximalSortSat <: AbstractSortSat end
struct DistalSortSat <: AbstractSortSat end

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
    #Sat(dist, db, convert(UInt32, root), P, C, cov)
    Sat(dist, db, convert(UInt32, root), C, cov)
end

@inline getpools(::Sat) = nothing
@inline database(sat::Sat) = sat.db
@inline database(sat::Sat, i) = sat.db[i]
@inline distance(sat::Sat) = sat.dist
@inline Base.length(sat::Sat) = length(sat.cov)

"""
    index!(
        sat::Sat;
        sortsat::AbstractSortSat=ProximalSortSat(),
        minleaf::Int=log2(ceil(database(sat)))
    )

Performs the indexing of the referenced dataset in the spatial access tree.

# Arguments
- `sat`: The metric data structure

# Keyword arguments
- `sortsat`: The strategy to create the spatial access tree, it heavily depends on the order of elements while it is build. It accepts:
   - `RandomSortSat()`: children are randomized.
   - `ProximalSortSat()`: classical approach, near elements are put first.
   - `DistalSortSat()`: recent approach, distant elements are put first.
- `minleaf`: Minimum number of children to perform a spatial access separation (half space partitioning)
"""
function index!(
        sat::Sat;
        sortsat::AbstractSortSat=ProximalSortSat(),
        minleaf::Int=log2(ceil(database(sat)))
    )
    n = length(sat)
    D = Vector{Tuple{Float32,UInt32}}(undef, n)
    p::UInt32 = sat.root
    sat.children[p] = collect(Iterators.flatten((1:p-1, p+1:n)))
    # sat.parents[p] = 0 # root has no parent
    sat.cov[p] = 0f0 # initializes cov for n = 1
    queue = UInt32[p]
    
    while length(queue) > 0
        p = pop!(queue)
        index_sat_neighbors!(sat, sat.children[p], D, p, sortsat, minleaf)
        for c in sat.children[p]
            C = sat.children[c]
            if C !== nothing
                push!(queue, c)
            end
        end
    end

    sat
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

function index_sat_neighbors!(sat::Sat, C::Nothing, D, p::UInt32, sortsat::AbstractSortSat)
    # do nothing
end

function index_sat_neighbors!(sat::Sat, C::Vector, D, p::UInt32, sortsat::AbstractSortSat, minleaf::Integer)
    # note: D is a cache of distances and objects, it is used in two ways in this function
    C = sat.children[p]
    n = length(C)

    resize!(D, n)
    parent = database(sat, p)
    dist = distance(sat)

    # computing distance to its parent (stored in D)
    sat.cov[p] = 0f0
    #L = Threads.SpinLock()

    for i in eachindex(C) #(i, c) in enumerate()
        c = C[i]
        d = evaluate(dist, parent, database(sat, c))
        D[i] = (convert(Float32, d), c)
        # not thread-safe:
        sat.cov[p] = max(sat.cov[p], d) # covering radius
    end

    T = typeof(sortsat)
    sort!(D, by=first)

    minleaf = min(n, minleaf)
    @inbounds for i in 1:minleaf
        (d_, i_) = D[i]
        D[i] = (-d_, i_)
    end

    if minleaf < n
        if T === RandomSortSat
            shuffle!(D)
        elseif T === DistalSortSat
            reverse!(D)
        end
    end

    # computing nearest neighbors of $child \in D$ (using previous D and storing the new set on D)
    empty!(C)

    for (d_, i_) in D
        d_ <= 0  && continue # negative distances encode mandatory leafs, see next outside-for-loop

        res = getknnresult(1)
        push!(res, p, d_)  # insert parent
        child = database(sat, i_)

        for j in C
            d = evaluate(dist, child, database(sat, j))
            push!(res, j, d)
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
    push!(res, p, dqp)

    if sat.children[p] !== nothing # inner node
        if length(res) < maxlength(res) || dqp < maximum(res) + sat.cov[p]
            for c in sat.children[p]
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


# this file is part of SpatialAccessTrees.jl

export Sat, RandomSortSat, ProximalSortSat, DistalSortSat, search, searchbatch, index!

abstract type AbstractSortSat end
struct RandomSortSat <: AbstractSortSat end
struct ProximalSortSat <: AbstractSortSat end
struct DistalSortSat <: AbstractSortSat end

struct Sat{DT<:SemiMetric,DBT<:AbstractDatabase} <: AbstractSearchIndex
    dist::DT
    db::DBT
    root::UInt32
    parents::Vector{UInt32}
    children::Vector{Union{Nothing,Vector{UInt32}}}
    cov::Vector{Float32} # leafs: d(parent, leaf), internal: max {d(parent, u) | u \in children(parent)}
end

function Sat(db::AbstractDatabase; dist::SemiMetric=L2Distance(), root=1)
    n = length(db)
    P = zeros(UInt32, n)
    C = Union{Nothing,Vector{UInt32}}[nothing for _ in 1:n]
    cov = Vector{Float32}(undef, n)
    Sat(dist, db, convert(UInt32, root), P, C, cov)
end

@inline getpools(::Sat) = nothing
@inline database(sat::Sat) = sat.db
@inline database(sat::Sat, i) = sat.db[i]
@inline distance(sat::Sat) = sat.dist
@inline Base.length(sat::Sat) = length(sat.cov)

function index!(sat::Sat; sortsat::AbstractSortSat=ProximalSortSat(), minleaf::Int=log2(ceil(database(sat))), randomize=false)
    n = length(sat)
    D = Vector{Tuple{Float32,UInt32}}(undef, n)
    p::UInt32 = sat.root
    sat.children[p] = collect(Iterators.flatten((1:p-1, p+1:n)))
    randomize && shuffle!(sat.children[p])
    sat.parents[p] = 0 # root has no parent
    sat.cov[p] = 0f0 # initializes cov for n = 1
    queue = UInt32[p]
    while length(queue) > 0
        p = pop!(queue)
        index_sat_neighbors!(sat, sat.children[p], D, p, sortsat)
        for c in sat.children[p]
            C = sat.children[c]
            if C !== nothing
                if length(C) >= minleaf
                    push!(queue, c)
                else
                    sat.cov[c] = abs(sum(sat.cov[i] for i in C))
                end
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

function index_sat_neighbors!(sat::Sat, C::Vector, D, p::UInt32, sortsat::AbstractSortSat)
    # note: D is a cache of distances and objects, it is used in two ways in this function
    n = length(sat.children[p])

    resize!(D, n)
    parent = database(sat, p)
    dist = distance(sat)

    # computing distance to its parent (stored in D)
    sat.cov[p] = -Inf32
    for (i, c) in enumerate(sat.children[p])
        d = evaluate(dist, parent, database(sat, c))
        D[i] = (convert(Float32, d), c)
        sat.cov[p] = max(sat.cov[p], d) # covering radius
    end

    T = typeof(sortsat)
    if T === RandomSortSat
        shuffle!(D)
    else
        sort!(D, by=first, rev=(T === DistalSortSat))
    end

    # computing nearest neighbors of $child \in D$ (using previous D and storing the new set on D)
    empty!(sat.children[p])
    #push!(C, last(D[1]))

    for (d_, i_) in D
        res = getknnresult(1)
        push!(res, p, d_)
        child = database(sat, i_)

        for j in sat.children[p]
            d = evaluate(dist, child, database(sat, j))
            push!(res, j, d)
        end

        nn = argmin(res)
        
        sat.cov[i_] = -minimum(res) # distance to parent, marked as negative
        sat.parents[i_] = nn

        if sat.children[nn] === nothing
            sat.children[nn] = UInt32[i_]
        else
            push!(sat.children[nn], i_)
        end
    end
end

isleaf(sat::Sat, i::Integer) = sat.cov[i] < 0
isinner(sat::Sat, i::Integer) = sat.cov[i] >= 0
isroot(sat::Sat, i::Integer) = sat.root == i

#=
function searchtree(sat::Sat, q, p::Integer, res::KnnResult)
    cost = 0
    r = sat.cov[p]

    if r >= 0  # inner node
        d = evaluate(distance(sat), q, database(sat, p))
        cost += 1
        push!(res, p, d)

        #if length(res) <= maxlength(res) #|| d < maximum(res) + r
            for c in sat.children[p]
                cost += searchtree(sat, q, c, res)
            end
        #end
    else
        #r = abs(r)
        #parent_ = sat.parents[p]
        #parentcov = convert(Float32, parent_ == 0 ? 0f0 : sat.cov[parent_])
        #if length(res) < maxlength(res) || r < parentcov + maximum(res)
            d = evaluate(distance(sat), q, database(sat, p))
            cost += 1
            push!(res, p, d)
        #end
    end

    cost
end
=#
function searchtree(sat::Sat, q, p::Integer, res::KnnResult)
    cost = 1
    d = evaluate(distance(sat), q, database(sat, p))
    push!(res, p, d)

    if sat.children[p] !== nothing # inner node
        if length(res) < maxlength(res) || d < maximum(res) + sat.cov[p]
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
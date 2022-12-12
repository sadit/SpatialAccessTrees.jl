# this file is part of SpatialAccessTrees.jl

export Sat, RandomSortSat, ProximalSortSat, DistalSortSat, search, searchbatch, index!

abstract type AbstractSortSat end
struct RandomSat <: AbstractSortSat end
struct ProximalSortSat <: AbstractSortSat end
struct DistalSortSat <: AbstractSortSat end

struct Sat{DT<:SemiMetric,DBT<:AbstractDatabase} <: AbstractSearchIndex
    dist::DT
    db::DBT
    parents::Vector{UInt32}
    children::Vector{Union{Nothing,Vector{UInt32}}}
    data::Vector{UInt32}
end

function Sat(db::AbstractDatabase; dist::SemiMetricL2Distance()=L2Distance())
    P = zeros(UInt32, length(db))
    C = Union{Nothing,Vector{UInt32}}[nothing for _ in eachindex(db)]
    D = Vector{UInt32}(undef, 0)
    sizehint!(D, length(db))
    Sat(dist, db, P, C, D)
end

@inline getpools(::Sat) = nothing
@inline database(sat::Sat) = sat.db
@inline database(sat::Sat, i) = sat.db[i]
@inline distance(sat::Sat) = sat.dist
@inline Base.length(sat::Sat) = length(sat.data)

function index!(sat::Sat, sortsat::AbstractSortSat=RandomSat())
    n = length(sat)
    T = NamedTuple{(:dist,:id), Tuple{Float32,UInt32}}
    D = Vector{T}(undef, n)
    p = 1
    sat.children[p] = collect(2:n)
    sat.parents[p] = 0
    _index_sat_neighbors!(sat, D, p, sortsat)
end

function index_sat_neighbors!(sat::Sat, D, p::Unsigned, sortsat::Union{ProximalSortSat,DistalSortSat})
    # note: D is a cache of distances and objects, it is used in two ways in this function
    C = sat.children[p]
    n = length(C)
    resize!(D, n)
    parent = database(sat, p)
    dist = distance(sat)

    # computing distance to its parent (stored in D)    
    for (i, c) in enumerate(C)
        d = evaluate(dist, parent, database(sat, c))
        D[i] = (d, c)
    end
    
    sort!(D, by=first, rev=(sortsat isa DistalSortSat))

    # computing nearest neighbors of $child \in D$ (using previous D and storing the new set on D)
    empty!(C)
    push!(C, D[1].id)
    for i in 2:n
        c = D[i]
        child = database(sat, c.id)
        res = getknnresult(1)
        push!(res, p, c.dist)
    
        for j in C
            d = evaluate(dist, child, database(sat, j))
            push!(res, j, d)
        end

        nn = argmin(res)
        if sat.children[nn] === nothing
            sat.children[nn] = UInt32[i]
        else
            push!(sat.children[nn], i)
        end
    end

    for c in C
        sat.parents[c] = p
        index_sat_neighbors!(sat, D, c, sortsat)
    end
end
# this file is part of SpatialAccessTrees.jl

export MultiSat, ParSat, allknn

struct ParSat{ST<:Sat} <: AbstractSearchIndex
    sat::ST
    parents::Vector{UInt32}
    depth::Int

    function ParSat(sat::S, depth::Int) where {S<:Sat}
        new{S}(sat, compute_parents(sat), depth)
    end
end

@inline getpools(::ParSat) = nothing
@inline database(psat::ParSat) = database(psat.sat)
@inline database(psat::ParSat, i) = database(psat.sat, i)
@inline distance(psat::ParSat) = distance(psat.sat)
@inline Base.length(psat::ParSat) = length(psat.sat)
search(psat::ParSat, q, res::KnnResult; pools=nothing) = search(psat.sat, q, res; pools)

function ith_par(parents, c, ith)
    j = 0
    while j < ith
        c = parents[c]
        j += 1
    end

    c
end

function compute_parents(sat::Sat)
    n = length(sat)
    P = Vector{UInt32}(undef, n)
    P[sat.root] = sat.root
    compute_parents(sat, P, sat.root)
    P
end

function compute_parents(sat::Sat, P::Vector, r::Integer)
    for c in sat.children[r]
        P[c] = r
        if sat.children[c] !== nothing
            compute_parents(sat, P, c)
        end
    end
end

function SimilaritySearch.allknn_single_search(psat::ParSat, i::Integer, k::Integer, pools)
    res = getknnresult(k, pools)
    q = database(psat, i)
    travelsat(psat.sat, q, ith_par(psat.parents, i, psat.depth), res)
    res
end

function travelsat(sat::Sat, q, p::Integer, res::KnnResult)
    dist = distance(sat)

    for c in sat.children[p]
        d = evaluate(dist, q, database(sat, c))
        push!(res, c, d)
        if sat.children[c] !== nothing && (length(res) < maxlength(res) || d < maximum(res) + sat.cov[c])
            travelsat(sat, q, c, res)
        end
    end
end

########## multisat

struct MultiSat{ST<:Sat} <: AbstractSearchIndex
    arr::Vector{ST}
    parents::Vector{Vector{UInt32}}
    depth::Int

    function MultiSat(arr, depth::Int)
        S = [sat for sat in arr]
        new{eltype(S)}(S, [compute_parents(sat) for sat in arr], depth)
    end
end

function Base.show(io::IO, msat::MultiSat)
    println(io, "{\n  ", typeof(msat), " size=", length(msat.arr), "\n")
    show(io, msat.arr)
    println(io, "}")
end

@inline getpools(::MultiSat) = nothing
@inline database(msat::MultiSat) = database(msat.arr[1])
@inline database(msat::MultiSat, i) = database(msat.arr[1], i)
@inline distance(msat::MultiSat) = distance(msat.arr[1])
@inline Base.length(msat::MultiSat) = length(msat.arr[1])

function SimilaritySearch.allknn_single_search(msat::MultiSat, i::Integer, k::Integer, pools)
    res = getknnresult(k, pools)
    vstate = reuse!(GlobalVisitedVertices[Threads.threadid()], length(msat.arr[1]))
    q = database(msat, i)

    for (sat, parents) in zip(msat.arr, msat.parents)
        if sat.children[i] === nothing
            travelmultisat(sat, q, ith_par(parents, i, msat.depth), res, vstate)
        else
            travelmultisat(sat, q, ith_par(parents, i, msat.depth), res, vstate)
        end
    end

    res
end

function travelmultisat(sat::Sat, q, p::Integer, res::KnnResult, vstate)
    dist = distance(sat)

    for c in sat.children[p]
        d = evaluate(dist, q, database(sat, c))
        !check_visited_and_visit!(vstate, convert(UInt64, p)) && push!(res, c, d)
        if sat.children[c] !== nothing && (length(res) < maxlength(res) || d < maximum(res) + sat.cov[c])
            travelmultisat(sat, q, c, res, vstate)
        end
    end
end

search(msat::MultiSat, q, res::KnnResult; pools=nothing) = search(msat.arr[1], q, res; pools)
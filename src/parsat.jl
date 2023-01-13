# this file is part of SpatialAccessTrees.jl

export ParSat, allknn

mutable struct ParSatConfig
    depth::Int32
    factor::Float32

    ParSatConfig(depth::Integer, factor::Real) = new(convert(Int32, depth), convert(Float32, factor))
end

struct ParSat{ST<:Sat} <: AbstractSearchIndex
    sat::ST
    parents::Vector{UInt32}
    config::ParSatConfig

    function ParSat(sat::S; depth::Integer=3, factor::Real=0.95) where {S<:Sat}
        new{S}(sat, compute_parents(sat), ParSatConfig(depth, factor))
    end
end

@inline getpools(::ParSat) = nothing
@inline database(psat::ParSat) = database(psat.sat)
@inline database(psat::ParSat, i) = database(psat.sat, i)
@inline distance(psat::ParSat) = distance(psat.sat)
@inline Base.length(psat::ParSat) = length(psat.sat)
search(psat::ParSat, q, res::KnnResult; pools=nothing) = search(psat.sat, q, res; pools)

"""
    ith_par(parents, c, ith)

Retrieves the ``i``-th parent of ``c``
"""
@inline function ith_par(parents, c, ith)
    j = 0
    @inbounds while j < ith
        c = parents[c]
        j += 1
    end

    c
end

"""
    compute_parents(sat::Sat) -> Vector{UInt32}

Computes the parents of all elements in the `sat` tree
"""
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

"""
    allknn_single_search(psat::ParSat, i::Integer, res::KnnResult, pools)

Solves a single query, ``i``th object, called from `allknn` function
"""
allknn_single_search(psat::ParSat, i::Integer, res::KnnResult, pools) = travelsat(psat, i, res, psat.config)

"""
    travelsat(psat::ParSat, i::Integer, res::KnnResult, c::ParSatConfig=psat.config)

Solves a single query, ``i``th object, specific solution for `allknn_single_search`, also used for
`runconfig`.
"""
function travelsat(psat::ParSat, i::Integer, res::KnnResult, c::ParSatConfig=psat.config)
    q = database(psat, i)
    p = ith_par(psat.parents, i, c.depth)
    cost = pruningsearchtree(psat.sat, q, p, res, c.factor)
    (; res, cost)
end

#=
function travelsat_(sat::Sat, q, p::Integer, res::KnnResult, factor::Float32)
    dist = distance(sat)
    cost = 0
    C = sat.children[p]

    @inbounds for c in C
        d = evaluate(dist, q, database(sat, c))
        push!(res, c, d)
        if sat.children[c] !== nothing && (length(res) < maxlength(res) || d < factor * maximum(res) + sat.cov[c])
            cost += travelsat_(sat, q, c, res, factor)
        end
    end

    cost + length(C)
end=#

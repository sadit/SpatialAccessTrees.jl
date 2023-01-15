# this file is part of SpatialAccessTrees.jl

export PrunParSat, allknn

mutable struct PrunParSatConfig
    depth::Int32
    factor::Float32

    PrunParSatConfig(depth::Integer, factor::Real) = new(convert(Int32, depth), convert(Float32, factor))
end

struct PrunParSat{ST<:Sat} <: AbstractSearchIndex
    sat::ST
    parents::Vector{UInt32}
    config::PrunParSatConfig

    function PrunParSat(sat::S; depth::Integer=3, factor::Real=0.95) where {S<:Sat}
        new{S}(sat, compute_parents(sat), PrunParSatConfig(depth, factor))
    end
end

@inline getpools(::PrunParSat) = nothing
@inline database(psat::PrunParSat) = database(psat.sat)
@inline database(psat::PrunParSat, i) = database(psat.sat, i)
@inline distance(psat::PrunParSat) = distance(psat.sat)
@inline Base.length(psat::PrunParSat) = length(psat.sat)
search(psat::PrunParSat, q, res::KnnResult; pools=nothing) = search(psat.sat, q, res; pools)

"""
    ith_par(parents, c, ith)

Retrieves the ``i``-th parent of ``c``
"""
@inline function ith_par(parents, c::Integer, ith::Integer)
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
    allknn_single_search(psat::PrunParSat, i::Integer, res::KnnResult, pools)

Solves a single query, ``i``th object, called from `allknn` function
"""
allknn_single_search(psat::PrunParSat, i::Integer, res::KnnResult, pools) = travelsat(psat, i, res, psat.config)

"""
    travelsat(psat::PrunParSat, i::Integer, res::KnnResult, c::PrunParSatConfig)

Solves a single query, ``i``th object, specific solution for `allknn_single_search`, also used for
`runconfig`.
"""
function travelsat(psat::PrunParSat, i::Integer, res::KnnResult, c::PrunParSatConfig)
    q = database(psat, i)
    p = ith_par(psat.parents, i, c.depth)
    cost = pruningsearchtree(psat.sat, q, p, res, c.factor)
    (; res, cost)
end

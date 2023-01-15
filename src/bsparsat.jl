# this file is part of SpatialAccessTrees.jl

export BeamSearchParSat, allknn

mutable struct BeamSearchParSatConfig
    depth::Int32
    bs::BeamSearch

    BeamSearchParSatConfig(depth::Integer, bs::BeamSearch) = new(convert(Int32, depth), bs)
end

struct BeamSearchParSat{ST<:Sat} <: AbstractSearchIndex
    sat::ST
    parents::Vector{UInt32}
    config::BeamSearchParSatConfig

    function BeamSearchParSat(sat::S; depth::Integer=3, bsize::Integer=4, Δ::Real=1.0) where {S<:Sat}
        c = BeamSearchParSatConfig(depth, BeamSearch(; bsize, Δ))
        new{S}(sat, compute_parents(sat), c)
    end
end

@inline getpools(::BeamSearchParSat) = nothing
@inline database(psat::BeamSearchParSat) = database(psat.sat)
@inline database(psat::BeamSearchParSat, i) = database(psat.sat, i)
@inline distance(psat::BeamSearchParSat) = distance(psat.sat)
@inline Base.length(psat::BeamSearchParSat) = length(psat.sat)
search(psat::BeamSearchParSat, q, res::KnnResult; pools=nothing) = search(psat.sat, q, res; pools)


"""
    allknn_single_search(psat::BeamSearchParSat, i::Integer, res::KnnResult, pools)

Solves a single query, ``i``th object, called from `allknn` function
"""
allknn_single_search(psat::BeamSearchParSat, i::Integer, res::KnnResult, pools) = beamsearch(psat, i, res, psat.config)
 
"""
    beamsearch(psat::BeamSearchParSat, i::Integer, res::KnnResult, c::BeamSearchParSatConfig)

Solves a single query, ``i``th object, specific solution for `allknn_single_search`, also used for
`runconfig`.
"""
function beamsearch(psat::BeamSearchParSat, i::Integer, res::KnnResult, c::BeamSearchParSatConfig)
    q = database(psat, i)
    p = ith_par(psat.parents, i, c.depth)
    beamsearch(psat.sat, p, c.bs.bsize, c.bs.Δ, q, res)
end

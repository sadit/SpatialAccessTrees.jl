# this file is part of SpatialAccessTrees.jl

export BeamSearchSat, optimize!

struct BeamSearchSat{ST<:Sat} <: AbstractSearchIndex
    bs::BeamSearch
    sat::ST
end

"""
    BeamSearchSat(sat::Sat; bsize=8, Δ=1f0)

Creates an approximate similarity sarch index based on sat and aggressive pruning; adapted from the paper
_A probabilistic spell for the curse of dimensionality_ (Chavez and Navarro, 2001).
It supports auto-tuning via  [`optimize!`](@ref).
"""
BeamSearchSat(sat::Sat; bsize=8, Δ=1f0) =
    BeamSearchSat(BeamSearch(; bsize, Δ), sat)

function Base.show(io::IO, bs::BeamSearchSat)
    println(io, "{\n  ", typeof(bs), " beam=", bs.bs, "\n")
    show(io, bs.sat)
    println(io, "}")
end

@inline getpools(::BeamSearchSat) = nothing
@inline database(bs::BeamSearchSat) = database(bs.sat)
@inline database(bs::BeamSearchSat, i) = database(bs.sat, i)
@inline distance(bs::BeamSearchSat) = distance(bs.sat)
@inline Base.length(bs::BeamSearchSat) = length(bs.sat)

function beamsearch(sat::Sat, root::UInt32, bsize::Int32, Δ::Float32, q, res::KnnResult)
    beam = reuse!(BeamKnnResult[Threads.threadid()], bsize)
    # beam = KnnResult(bs.bsize)
    dist = distance(sat)
    cost = 1
    d = evaluate(dist, q, database(sat, root))
    push_item!(res, IdWeight(root, d))
    sat.children[root] !== nothing && push_item!(beam, IdWeight(root, d))

    # bsize = maxlength(beam)
    sp = 1

    @inbounds while sp <= length(beam)
        i_ = beam[sp].id
        sp += 1
        C = sat.children[i_]::Vector{UInt32}
        for c in C
            d = evaluate(dist, q, database(sat, c))
            cost += 1
            push_item!(res, IdWeight(c, d))
        
            if sat.children[c] !== nothing && d <= Δ * maximum(res)
                push_item!(beam, IdWeight(c, d), sp)
            end
        end
    end

    (; res, cost)
end

search(bs::BeamSearchSat, q, res::KnnResult; pools=nothing) = beamsearch(bs.sat, bs.sat.root, bs.bs.bsize, bs.bs.Δ, q, res)

## Optimization

optimization_space(::BeamSearchSat) =
    BeamSearchSpace(;
        bsize = 8:16:64,
        Δ = [0.8, 1.0, 1.3, 1.5],
        bsize_scale = (s=1.5, p1=0.5, p2=0.5, lower=4, upper=256),
        Δ_scale = (s=1.2, p1=0.5, p2=0.5, lower=0.5, upper=3.0)
    )

function optimize!(
        index::BeamSearchSat,
        kind::ErrorFunction=MinRecall(0.9),
        space::AbstractSolutionSpace=optimization_space(index);
        queries=nothing,
        ksearch=10,
        numqueries=64,
        initialpopulation=16,
        minbatch=0,
        maxpopulation=16,
        bsize=8,
        mutbsize=16,
        crossbsize=8,
        tol=-1.0,
        maxiters=16,
        verbose=false,
    )

    SimilaritySearch.optimize_index!(index, kind, space; queries, ksearch, numqueries, initialpopulation, minbatch,
        maxpopulation, bsize, mutbsize, crossbsize, tol, maxiters, verbose)
end

function setconfig!(bs::BeamSearch, index::BeamSearchSat, perf)
    index.bs.bsize = bs.bsize
    index.bs.Δ = bs.Δ
end

function runconfig(bs::BeamSearch, index::BeamSearchSat, q, res::KnnResult, pools)
    beamsearch(index.sat, index.sat.root, bs.bsize, bs.Δ, q, res)
end
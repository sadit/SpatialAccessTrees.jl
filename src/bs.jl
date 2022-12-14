# this file is part of SpatialAccessTrees.jl

export BeamSearchSat, optimize!

struct BeamSearchSat{ST<:Sat} <: AbstractSearchIndex
    bs::BeamSearch
    sat::ST
end

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

function beamsearch(sat::Sat, bsize::Int32, Δ::Float32, q, res::KnnResult)
    beam = reuse!(BeamKnnResult[Threads.threadid()], bsize)
    # beam = KnnResult(bs.bsize)
    dist = distance(sat)
    root = sat.root
    cost = 1
    d = evaluate(dist, q, database(sat, root))
    push!(res, root, d)
    sat.children[root] !== nothing && push!(beam, root, d)

    
    while length(beam) > 0
        (i_, _) = popfirst!(beam)
        for c in sat.children[i_]
            d = evaluate(dist, q, database(sat, c))
            cost += 1
            push!(res, c, d)
        
            if sat.children[c] !== nothing && d <= Δ * maximum(res)
                push!(beam, c, d)
            end
        end
    end

    (; res, cost)
end

search(bs::BeamSearchSat, q, res::KnnResult; pools=nothing) = beamsearch(bs.sat, bs.bs.bsize, bs.bs.Δ, q, res)

## Optimization

optimization_space(::BeamSearchSat) =
    BeamSearchSpace(;
        Δ = [0.8, 1.0, 1.3, 1.5],
        bsize = 8:16:64,
        bsize_scale = (s=1.5, p1=0.5, p2=0.5, lower=4, upper=256),
        Δ_scale = (s=1.1, p1=0.5, p2=0.5, lower=0.5, upper=3.0)
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
    beamsearch(index.sat, bs.bsize, bs.Δ, q, res)
end
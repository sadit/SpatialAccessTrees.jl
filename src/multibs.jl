# this file is part of SpatialAccessTrees.jl

export BeamSearchMultiSat

struct BeamSearchMultiSat{ST<:Sat} <: AbstractSearchIndex
    bs::BeamSearch
    sat::Vector{ST}
end

BeamSearchMultiSat(sat::Vector{SatType}; bsize=8, Δ=1f0) where {SatType<:Sat} =
    BeamSearchMultiSat(BeamSearch(; bsize, Δ), sat)

function Base.show(io::IO, bs::BeamSearchMultiSat)
    println(io, "{\n  ", typeof(bs), " beam=", bs.bs, " size=", length(bs.sat), "\n")
    show(io, bs.sat)
    println(io, "}")
end

@inline getpools(::BeamSearchMultiSat) = nothing
@inline database(bs::BeamSearchMultiSat) = database(bs.sat[1])
@inline database(bs::BeamSearchMultiSat, i) = database(bs.sat[1], i)
@inline distance(bs::BeamSearchMultiSat) = distance(bs.sat[1])
@inline Base.length(bs::BeamSearchMultiSat) = length(bs.sat[1])

function beamsearch(satarr::Vector{SatType}, bsize::Int32, Δ::Float32, q, res::KnnResult) where {SatType<:Sat}
    beam = reuse!(BeamKnnResult[Threads.threadid()], bsize)
    sat = satarr[1]
    vstate = reuse!(GlobalVisitedVertices[Threads.threadid()], length(sat))
    dist = distance(sat)
    root = sat.root
    cost = 1
    d = evaluate(dist, q, database(sat, root))
    push!(res, root, d)
    sat.children[root] !== nothing && push!(beam, root, d)
    visit!(vstate, convert(UInt64, root))

    round_robin, m = 0, length(satarr)
    bsize = maxlength(beam)
    sp = 1

    @inbounds while sp <= length(beam)
        i_ = getid(beam, sp)
        sp += 1
        round_robin = ifelse(round_robin == m, 1, round_robin + 1)
        # round_robin = (round_robin > m) ? 1 : round_robin + 1
        sat = satarr[round_robin]

        C = sat.children[i_]
        C === nothing && continue
        for c in C
            d = evaluate(dist, q, database(sat, c))
            cost += 1
            push!(res, c, d)
        
            if sat.children[c] !== nothing && d <= Δ * maximum(res)
                # push!(beam, c, d)
                check_visited_and_visit!(vstate, convert(UInt64, c)) && continue
                push!(beam, c, d; sp, k=bsize+sp)
            end
        end
    end

    (; res, cost)
end

search(bs::BeamSearchMultiSat, q, res::KnnResult; pools=nothing) = beamsearch(bs.sat, bs.bs.bsize, bs.bs.Δ, q, res)

## Optimization

optimization_space(::BeamSearchMultiSat) =
    BeamSearchSpace(;
        Δ = [0.8, 1.0, 1.3, 1.5],
        bsize = 8:16:64,
        bsize_scale = (s=1.5, p1=0.5, p2=0.5, lower=4, upper=256),
        Δ_scale = (s=1.1, p1=0.5, p2=0.5, lower=0.5, upper=3.0)
    )

function optimize!(
        index::BeamSearchMultiSat,
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

function setconfig!(bs::BeamSearch, index::BeamSearchMultiSat, perf)
    index.bs.bsize = bs.bsize
    index.bs.Δ = bs.Δ
end

function runconfig(bs::BeamSearch, index::BeamSearchMultiSat, q, res::KnnResult, pools)
    beamsearch(index.sat, bs.bsize, bs.Δ, q, res)
end
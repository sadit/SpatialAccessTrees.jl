# this file is part of SpatialAccessTrees.jl
export BeamSearchSat

struct BeamSearchSat{ST<:Sat} <: AbstractSearchIndex
    sat::ST
    bsize::UInt32
    Δ::Float32
end

BeamSearchSat(sat::Sat; bsize=8, Δ=1f0) =
    BeamSearchSat(sat, convert(UInt32, bsize), convert(Float32, Δ))

function Base.show(io::IO, bs::BeamSearchSat)
    println(io, "{\n  ", typeof(bs), " bsize=", length(bs.bsize), "\n")
    show(io, bs)
    println(io, "}")
end

@inline getpools(::BeamSearchSat) = nothing
@inline database(bs::BeamSearchSat) = database(bs.sat)
@inline database(bs::BeamSearchSat, i) = database(bs.sat, i)
@inline distance(bs::BeamSearchSat) = distance(bs.sat)
@inline Base.length(bs::BeamSearchSat) = length(bs.sat)

function beamsearch(sat::Sat, bsize::UInt32, Δ::Float32, q, res::KnnResult)
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

search(bs::BeamSearchSat, q, res::KnnResult; pools=nothing) = beamsearch(bs.sat, bs.bsize, bs.Δ, q, res)

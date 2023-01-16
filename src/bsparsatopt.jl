# this file is part of SpatialAccessTrees.jl

export BeamSearchParSatConfigSpace

@with_kw struct BeamSearchParSatConfigSpace <: AbstractSolutionSpace
    depth = 1:2:7
    depth_scale = (s=1.5, p1=0.5, p2=0.5, lower=1, upper=15)
    bs = BeamSearchSpace(;
        bsize = 8:16:64,
        Δ = [0.8, 1.0, 1.3, 1.5],
        bsize_scale = (s=1.5, p1=0.5, p2=0.5, lower=4, upper=256),
        Δ_scale = (s=1.2, p1=0.5, p2=0.5, lower=0.5, upper=3.0)
    )
end

Base.hash(c::BeamSearchParSatConfig) = hash((c.depth, hash(c.bs)))
Base.isequal(a::BeamSearchParSatConfig, b::BeamSearchParSatConfig) = (a.depth == b.depth) && (a.bs == b.bs)
Base.eltype(::BeamSearchParSatConfigSpace) = BeamSearchParSatConfig
Base.rand(space::BeamSearchParSatConfigSpace) = BeamSearchParSatConfig(rand(space.depth), rand(space.bs))

function combine(a::BeamSearchParSatConfig, b::BeamSearchParSatConfig)
    depth = ceil(Int, (a.depth + b.depth) / 2)
    bs = combine(a.bs, b.bs)
    BeamSearchParSatConfig(depth, bs)
end

function mutate(space::BeamSearchParSatConfigSpace, c::BeamSearchParSatConfig, iter)
    depth = SearchModels.scale(c.depth; space.depth_scale...)
    bs = mutate(space.bs, c.bs, iter)
    BeamSearchParSatConfig(depth, bs)
end

optimization_space(psat::BeamSearchParSat) = BeamSearchParSatConfigSpace()

function setconfig!(config::BeamSearchParSatConfig, psat::BeamSearchParSat, perf)
    psat.config.depth = config.depth
    psat.config.bs = config.bs
end

function runconfig0(config::BeamSearchParSatConfig, psat::BeamSearchParSat, queries::SubDatabase, i::Integer, res::KnnResult, pools)
    beamsearch(psat, convert(UInt32, queries.map[i]), res, config)
end

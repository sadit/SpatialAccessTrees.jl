# this file is part of SpatialAccessTrees.jl

export PrunParSatConfigSpace

@with_kw struct PrunParSatConfigSpace <: AbstractSolutionSpace
    depth = 1:2:7
    factor = 0.3:0.1:0.99
    depth_scale = (s=1.5, p1=0.5, p2=0.5, lower=1, upper=15)
    factor_scale = (s=1.07, p1=0.5, p2=0.5, lower=0.1, upper=0.9999)  # that should work in most datasets
end

Base.hash(c::PrunParSatConfig) = hash((c.depth, round(c.factor, digits=2)))
Base.isequal(a::PrunParSatConfig, b::PrunParSatConfig) = (a.depth == b.depth) && (a.factor == b.factor)
Base.eltype(::PrunParSatConfigSpace) = PrunParSatConfig
Base.rand(space::PrunParSatConfigSpace) = PrunParSatConfig(rand(space.depth), rand(space.factor))

function combine(a::PrunParSatConfig, b::PrunParSatConfig)
    depth = ceil(Int, (a.depth + b.depth) / 2)
    factor = (a.factor + b.factor) / 2
    PrunParSatConfig(depth, factor)
end

function mutate(space::PrunParSatConfigSpace, c::PrunParSatConfig, iter)
    depth = SearchModels.scale(c.depth; space.depth_scale...)
    factor = SearchModels.scale(c.factor; space.factor_scale...)
    PrunParSatConfig(depth, factor)
end

optimization_space(index::PrunParSat) = PrunParSatConfigSpace()

function setconfig!(config::PrunParSatConfig, psat::PrunParSat, perf)
    psat.config.depth = config.depth
    psat.config.factor = config.factor
end

function runconfig0(config::PrunParSatConfig, psat::PrunParSat, queries::SubDatabase, i::Integer, res::KnnResult, pools)
    travelsat(psat, queries.map[i], res, config)
end

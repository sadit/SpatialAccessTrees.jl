# this file is part of SpatialAccessTrees.jl

export ParSatConfigSpace

@with_kw struct ParSatConfigSpace <: AbstractSolutionSpace
    depth = 1:2:7
    factor = 0.3:0.03:0.99
    depth_scale = (s=1.5, p1=0.5, p2=0.5, lower=1, upper=15)
    factor_scale = (s=1.07, p1=0.5, p2=0.5, lower=0.1, upper=0.9999)  # that should work in most datasets
end

Base.hash(c::ParSatConfig) = hash((c.depth, round(c.factor, digits=4)))
Base.isequal(a::ParSatConfig, b::ParSatConfig) = (a.depth == b.depth) && (a.factor == b.factor)
Base.eltype(::ParSatConfigSpace) = ParSatConfig
Base.rand(space::ParSatConfigSpace) = ParSatConfig(rand(space.depth), rand(space.factor))

function combine(a::ParSatConfig, b::ParSatConfig)
    depth = ceil(Int, (a.depth + b.depth) / 2)
    factor = (a.factor + b.factor) / 2
    ParSatConfig(depth, factor)
end

function mutate(space::ParSatConfigSpace, c::ParSatConfig, iter)
    depth = SearchModels.scale(c.depth; space.depth_scale...)
    factor = SearchModels.scale(c.factor; space.factor_scale...)
    ParSatConfig(depth, factor)
end

optimization_space(index::ParSat) = ParSatConfigSpace()

function setconfig!(config::ParSatConfig, psat::ParSat, perf)
    psat.config.depth = config.depth
    psat.config.factor = config.factor
end

function runconfig0(config::ParSatConfig, psat::ParSat, queries::SubDatabase, i::Integer, res::KnnResult, pools)
    travelsat(psat, queries.map[i], res, config)
end

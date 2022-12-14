# this file is part of SpatialAccessTrees.jl

export PruningSatFactorSpace

@with_kw struct PruningSatFactorSpace <: AbstractSolutionSpace
    factor = 0.2:0.1:0.9
    factor_scale = (s=1.07, p1=0.5, p2=0.5, lower=0.1, upper=0.9999)  # that should work in most datasets
end

Base.hash(c::PruningSatFactor) = hash(round(c.factor, digits=4))
Base.isequal(a::PruningSatFactor, b::PruningSatFactor) = a.factor == b.factor
Base.eltype(::PruningSatFactorSpace) = PruningSatFactor
Base.rand(space::PruningSatFactorSpace) = PruningSatFactor(rand(space.factor))

function combine(a::PruningSatFactor, b::PruningSatFactor)
    x = ceil(Int, (a.factor + b.factor) / 2)
    PruningSatFactor(x)
end

function mutate(space::PruningSatFactorSpace, c::PruningSatFactor, iter)
    x = SearchModels.scale(c.factor; space.factor_scale...)
    PruningSatFactor(x)
end

optimization_space(index::PruningSat) = PruningSatFactorSpace()

function setconfig!(config::PruningSatFactor, index::PruningSat, perf)
    index.pruning.factor = config.factor
end

function runconfig(config::PruningSatFactor, index::PruningSat, q, res::KnnResult, pools)
    cost = pruningsearchtree(index.sat, q, index.sat.root, res, config.factor)
    (; res, cost)
end

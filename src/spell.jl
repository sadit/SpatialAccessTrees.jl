# this file is part of SpatialAccessTrees.jl

export PruningSat

mutable struct PruningSatFactor
    factor::Float32
end

struct PruningSat{ST<:Sat} <: AbstractSearchIndex
    pruning::PruningSatFactor
    sat::ST
end

"""
    PruningSat(sat::Sat; factor=0.9)

Creates an approximate similarity sarch index based on sat and aggressive pruning; adapted from the paper
_A probabilistic spell for the curse of dimensionality_ (Chavez and Navarro, 2001).
It supports auto-tuning via  [`optimize!`](@ref).
"""
PruningSat(sat::Sat; factor=0.9) = PruningSat(PruningSatFactor(convert(Float32, factor)), sat)

function Base.show(io::IO, psat::PruningSat)
    println(io, "{\n  ", typeof(psat), " pruning=", psat.pruning, "\n")
    show(io, psat.sat)
    println(io, "}")
end

@inline getpools(::PruningSat) = nothing
@inline database(psat::PruningSat) = database(psat.sat)
@inline database(psat::PruningSat, i) = database(psat.sat, i)
@inline distance(psat::PruningSat) = distance(psat.sat)
@inline Base.length(psat::PruningSat) = length(psat.sat)


explore_node!(sat::Sat, q, ::Nothing, res::KnnResult, queue::Vector) = 0

function explore_node!(sat::Sat, q, C::Vector, res::KnnResult, queue::Vector)
    cost = 0
    dist = distance(sat)
    for c in C
        if sat.children[c] === nothing
            d = evaluate(dist, q, database(sat, c))
            cost += 1
            push_item!(res, IdWeight(c, d))
        else
            push!(queue, c)
        end
    end

    cost
end

function pruningsearchtree(sat::Sat, q, p::Integer, res::KnnResult, factor::Float32)
    cost = 0
    queue = getsearchqueue()
    push!(queue, p)
    dist = distance(sat)
    
    @inbounds while length(queue) > 0
        p = pop!(queue)
        dqp = evaluate(dist, q, database(sat, p))
        cost += 1
        push_item!(res, IdWeight(p, dqp))

        if length(res) < maxlength(res) || dqp < factor * maximum(res) + sat.cov[p]
            cost += explore_node!(sat, q, sat.children[p], res, queue)
        end
    end

    cost
end

#=function pruningsearchtree(sat::Sat, q, p::Integer, res::KnnResult, factor::Float32)
    cost = 1
    dist = distance(sat)
    dqp = evaluate(dist, q, database(sat, p))
    push_item!(res, IdWeight(p, dqp))

    if sat.children[p] !== nothing # inner node
        if length(res) < maxlength(res) || dqp < factor * maximum(res) + sat.cov[p]
            for c in sat.children[p]
                cost += pruningsearchtree(sat, q, c, res, factor)
            end
        end
    end

    cost
end=#

function search(psat::PruningSat, q::T, res::KnnResult; pools=nothing) where T
    cost = pruningsearchtree(psat.sat, q, psat.sat.root, res, psat.pruning.factor)
    (; res, cost)
end


#- `pruning_factor`: factor for aggressive pruning candidates, valid range: ``0 < pruning_factor \\leq 1``, an exact search is made when `pruning_factor=1`.
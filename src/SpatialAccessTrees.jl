# this file is part of SpatialAccessTrees.jl

module SpatialAccessTrees

using SimilaritySearch, Parameters, SearchModels, Random, Polyester

import SimilaritySearch:
    search, index!, getpools, getknnresult, database, distance,
    optimization_space, optimize!, setconfig!, runconfig, mutate, combine

using SimilaritySearch:
    ErrorFunction, AbstractSolutionSpace, GlobalVisitedVertices, visit!, check_visited_and_visit!, getminbatch


const BeamKnnResult = [KnnResult(32)]  # see __init__ function
const SearchQueue = [UInt32[]]

function __init__()
    for _ in 2:Threads.nthreads()
        push!(BeamKnnResult, KnnResult(32))
        push!(SearchQueue, copy(first(SearchQueue)))
    end
end

include("sat.jl")
include("bs.jl")
include("spell.jl")
include("spellopt.jl")
include("multibs.jl")
include("permsat.jl")

end
 
# this file is part of SpatialAccessTrees.jl

module SpatialAccessTrees

using SimilaritySearch, Parameters, SearchModels, Random, Polyester

import SimilaritySearch:
    search, index!, getpools, getknnresult, database, distance,
    optimization_space, optimize!, allknn_single_search, setconfig!, runconfig0, runconfig, mutate, combine

using SimilaritySearch:
    ErrorFunction, AbstractSolutionSpace, GlobalVisitedVertices, visit!, check_visited_and_visit!, getminbatch


const BeamKnnResult = [KnnResult(32)]  # see __init__ function
const SearchQueue = [Vector{UInt32}(undef, 64)]

function __init__()
    for _ in 2:Threads.nthreads()
        push!(BeamKnnResult, KnnResult(32))
        push!(SearchQueue, copy(first(SearchQueue)))
    end
end

function getsearchqueue()
    queue = SearchQueue[Threads.threadid()]
    empty!(queue)
    queue
end

include("sat.jl")
include("bs.jl")
include("spell.jl")
include("spellopt.jl")
include("multibs.jl")
include("permsat.jl")
include("prunparsat.jl")
include("prunparsatopt.jl")
include("bsparsat.jl")
include("bsparsatopt.jl")

end
 
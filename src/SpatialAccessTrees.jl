module SpatialAccessTrees

using SimilaritySearch, Parameters, SearchModels, Random

import SimilaritySearch:
    search, index!, getpools, getknnresult, database, distance,
    optimization_space, optimize!, setconfig!, runconfig, mutate, combine

using SimilaritySearch:
    ErrorFunction, AbstractSolutionSpace


const BeamKnnResult = [KnnResult(32)]  # see __init__ function

function __init__()
    for _ in 2:Threads.nthreads()
        push!(BeamKnnResult, KnnResult(32))
    end
end

include("sat.jl")
include("bs.jl")
include("spell.jl")
include("spellopt.jl")

end
 
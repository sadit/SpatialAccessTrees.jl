module SpatialAccessTrees

using SimilaritySearch, Random
import SimilaritySearch: search, index!, evaluate, getpools, getknnresult, database, distance

const BeamKnnResult = [KnnResult(32)]  # see __init__ function

function __init__()
    for _ in 2:Threads.nthreads()
        push!(BeamKnnResult, KnnResult(32))
    end
end

include("sat.jl")
include("bs.jl")
end
 
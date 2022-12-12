using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees.jl" begin
    dim = 4
    db = MatrixDatabase(rand(Float32, dim, 10_000))
    queries = MatrixDatabase(rand(Float32, dim, 100))

    sat = SpatialAccessTrees(dist, db)
    index!(sat)
end

using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees.jl" begin
    dim = 2
    dist = L2Distance()
    n = 1000
    m = 100
    k = 4
    db = MatrixDatabase(rand(Float32, dim, n))
    queries = MatrixDatabase(rand(Float32, dim, m))
    E = ExhaustiveSearch(; dist, db)
    goldsearchtime = @elapsed Igold, Dgold = searchbatch(E, queries, k)

    for sortsat in [ProximalSortSat(), DistalSortSat(), RandomSortSat()]
        minleaf=4
        sat = Sat(db; dist)
        @show index!(sat; sortsat, minleaf)
        # checking that the database size and the number of inserted elements is consistent
        @test n == 1 + sum(length(C) for C in sat.children if C !== nothing)
        # checking that cov is consistent with children
        @test all((C === nothing ? sat.cov[i] < 0 : sat.cov[i] > 0) for (i, C) in enumerate(sat.children))
        searchtime = @elapsed Isat, Dsat = searchbatch(sat, queries, k)
        @show "---------" sortsat macrorecall(Igold, Isat) searchtime goldsearchtime
    end
end

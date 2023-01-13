# this file is part of SpatialAccessTrees.jl

using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees" begin
    n, k, dim, minleaf, dist = 10^4, 15, 8, 2, L2Distance()
    db = StrideMatrixDatabase(rand(Float32, dim, n))
    E = ExhaustiveSearch(; dist, db)
    bruteforcesearchtime = @elapsed Igold, _ = allknn(E, k)
    numqueries = 16  # numqueries for optimize!
    verbose = false
    sortsat = RandomSortSat()

    let
        buildtime = @elapsed psat = ParSat(index!(Sat(db; dist); sortsat, minleaf))
        blist = optimize!(psat, MinRecall(0.9); verbose=false, numqueries)
        println(stderr, first(blist))
        #buildtime = @elapsed psat = index!(Sat(db; dist); sortsat, minleaf)
        I, _ = allknn(psat, k)
        searchtime = @elapsed I, _ = allknn(psat, k)
        @show bruteforcesearchtime bruteforcesearchtime / searchtime
        @show searchtime buildtime macrorecall(Igold, I)
    end
end
# this file is part of SpatialAccessTrees.jl

using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees" begin
    n, k, dim, minleaf, dist = 10^5, 15, 8, 2, L2Distance()
    db = StrideMatrixDatabase(rand(Float32, dim, n))
    E = ExhaustiveSearch(; dist, db)
    bruteforcesearchtime = @elapsed Igold, _ = allknn(E, k)
    numqueries = 16  # numqueries for optimize!
    verbose = false
    sortsat = RandomSortSat()

    let
        buildtime = @elapsed psat = ParSat(index!(Sat(db; dist); sortsat, minleaf), 4)
        #buildtime = @elapsed psat = index!(Sat(db; dist); sortsat, minleaf)
        I, _ = allknn(psat, k)
        searchtime = @elapsed I, _ = allknn(psat, k)
        @show bruteforcesearchtime bruteforcesearchtime / searchtime
        @show searchtime buildtime macrorecall(Igold, I)
    end

   #= let
        buildtime = @elapsed msat = MultiSat([index!(Sat(db; dist, root=rand(1:n)); sortsat, minleaf) for _ in 1:5], 5)
        I, _ = allknn(msat, k)
        searchtime = @elapsed I, _ = allknn(msat, k)
        @show bruteforcesearchtime bruteforcesearchtime / searchtime
        @show searchtime buildtime macrorecall(Igold, I)
    end=#
end
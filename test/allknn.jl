# this file is part of SpatialAccessTrees.jl

using SpatialAccessTrees, SimilaritySearch
using Test
#=
using JET
let
    n, k, dim, minleaf, dist = 10^4, 15, 8, 2, L2Distance()
    db = StrideMatrixDatabase(rand(Float32, dim, n))
    E = ExhaustiveSearch(; dist, db)
    bruteforcesearchtime = @elapsed Igold, _ = allknn(E, k)
    numqueries = 16  # numqueries for optimize!
    verbose = false
    sortsat = RandomSortSat()

## let
##     buildtime = @elapsed psat = PrunParSat(index!(Sat(db; dist); sortsat, minleaf))
##     blist = optimize!(psat, MinRecall(0.9); verbose=false, numqueries)
##     println(stderr, first(blist))
##     @report_opt allknn(psat, k)
## end

let
        buildtime = @elapsed psat = BeamSearchParSat(index!(Sat(db; dist); sortsat, minleaf))
        blist = optimize!(psat, MinRecall(0.9); verbose=false, numqueries)
        println(stderr, first(blist))
        #buildtime = @elapsed psat = index!(Sat(db; dist); sortsat, minleaf)
        @report_opt allknn(psat, k)
end
end
=#

@testset "SpatialAccessTrees" begin
    n, k, dim, minleaf, dist = 10^4, 15, 8, 2, L2Distance()
    db = StrideMatrixDatabase(rand(Float32, dim, n))
    E = ExhaustiveSearch(; dist, db)
    bruteforcesearchtime = @elapsed Igold, _ = allknn(E, k)
    numqueries = 16  # numqueries for optimize!
    verbose = false
    sortsat = RandomSortSat()
    println(stderr, "================= PrunParSat ================")
    let
        buildtime = @elapsed psat = PrunParSat(index!(Sat(db; dist); sortsat, minleaf))
        blist = optimize!(psat, MinRecall(0.9); verbose=false, numqueries)
        println(stderr, first(blist))
        #buildtime = @elapsed psat = index!(Sat(db; dist); sortsat, minleaf)
        I, _ = allknn(psat, k)
        searchtime = @elapsed I, _ = allknn(psat, k)
        @show bruteforcesearchtime bruteforcesearchtime / searchtime
        @show searchtime buildtime macrorecall(Igold, I)
    end

    println(stderr, "================= BeamSearchParSat ================")
    let
        buildtime = @elapsed psat = BeamSearchParSat(index!(Sat(db; dist); sortsat, minleaf))
        blist = optimize!(psat, MinRecall(0.9); verbose=false, numqueries)
        println(stderr, first(blist))
        #buildtime = @elapsed psat = index!(Sat(db; dist); sortsat, minleaf)
        I, _ = allknn(psat, k)
        searchtime = @elapsed I, _ = allknn(psat, k)
        @show bruteforcesearchtime bruteforcesearchtime / searchtime
        @show searchtime buildtime macrorecall(Igold, I)
    end
end
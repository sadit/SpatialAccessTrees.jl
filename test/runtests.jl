using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees" begin
    dim = 8
    dist = L2Distance()
    n = 10_000
    m = 64
    k = 4
    db = MatrixDatabase(rand(Float32, dim, n))
    queries = MatrixDatabase(rand(Float32, dim, m))
    E = ExhaustiveSearch(; dist, db)
    goldsearchtime = @elapsed Igold, Dgold = searchbatch(E, queries, k)
    minleaf=4

    for sortsat in [ProximalSortSat(), DistalSortSat(), RandomSortSat()]
        @info "============================= $sortsat minleaf=$minleaf ================"
        
        sat = Sat(db; dist)
        @time index!(sat; sortsat, minleaf)
        @show sat
        # checking that the database size and the number of inserted elements is consistent
        @test n == 1 + sum(length(C) for C in sat.children if C !== nothing)
        # checking that cov is consistent with children
        @test all((C === nothing ? sat.cov[i] < 0 : sat.cov[i] > 0) for (i, C) in enumerate(sat.children))
        searchtime = @elapsed Isat, Dsat = searchbatch(sat, queries, k)
        recall = macrorecall(Igold, Isat)
        @test recall >= 0.9999

        asat = Sat(sat; pruning_factor=0.4)
        asearchtime = @elapsed Ia, Da = searchbatch(asat, queries, k)
        arecall = macrorecall(Igold, Ia)

        bsat = BeamSearchSat(sat)
        optimize!(bsat, MinRecall(0.5), verbose=true)
        bsearchtime = @elapsed Ib, Db = searchbatch(bsat, queries, k)
        @info bsat
        @time Ib, Db = searchbatch(bsat, queries, k)
        brecall = macrorecall(Igold, Ib)

        @show "------------ "
        @show " exact:" recall (goldsearchtime / searchtime)
        @show " probabilistic spell:" arecall (goldsearchtime / asearchtime)
        @show " beam search:" brecall (goldsearchtime / bsearchtime)
        break
    end
end

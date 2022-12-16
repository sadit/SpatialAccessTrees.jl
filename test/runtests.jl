using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees" begin
    dim = 8
    dist = L2Distance()
    n = 100_000
    m = 64
    k = 4
    db = MatrixDatabase(rand(Float32, dim, n))
    queries = MatrixDatabase(rand(Float32, dim, m))
    E = ExhaustiveSearch(; dist, db)
    bruteforcesearchtime = @elapsed Igold, Dgold = searchbatch(E, queries, k)
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
        Isat, Dsat = searchbatch(sat, queries, k)
        searchtime = @elapsed Isat, Dsat = searchbatch(sat, queries, k)
        recall = macrorecall(Igold, Isat)
        @test recall >= 0.9999

        asat = PruningSat(sat)
        optimize!(asat, MinRecall(0.9))
        Ia, _ = searchbatch(asat, queries, k)
        asearchtime = @elapsed Ia, _ = searchbatch(asat, queries, k)
        arecall = macrorecall(Igold, Ia)
        @test arecall > 0.6 ## it should be close to 0.9, but anyway errors will show very low recall

        bsat = BeamSearchSat(sat)
        optimize!(bsat, MinRecall(0.9), verbose=false)
        @show bsat.bs
        Ib, _ = searchbatch(bsat, queries, k)
        bsearchtime = @elapsed Ib, _ = searchbatch(bsat, queries, k)
        brecall = macrorecall(Igold, Ib)
        @test brecall > 0.6

        csat = BeamSearchMultiSat([index!(Sat(db; dist, root=rand(1:n)); sortsat, minleaf) for _ in 1:8])
        # @show [s.root for s in csat.sat]
        @show optimize!(csat, MinRecall(0.9), verbose=false)
        @show csat.bs
        Ic, _ = searchbatch(csat, queries, k)
        csearchtime = @elapsed Ic, _ = searchbatch(csat, queries, k)
        crecall = macrorecall(Igold, Ic)
        # @test crecall > 0.6

        @info "------------ brute force: $bruteforcesearchtime "
        @info " exact:" recall searchtime (bruteforcesearchtime / searchtime)
        @info " probabilistic spell:" arecall asearchtime (bruteforcesearchtime / asearchtime)
        @info " beam search:" brecall bsearchtime (bruteforcesearchtime / bsearchtime)
        @info " beam search - multi sat:" crecall csearchtime (bruteforcesearchtime / csearchtime)
        break
    end
end

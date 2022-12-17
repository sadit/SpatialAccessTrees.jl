using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees" begin
    dim = 4
    dist = L2Distance()
    n = 100_000
    m = 100
    k = 10
    db = StrideMatrixDatabase(rand(Float32, dim, n))
    queries = StrideMatrixDatabase(rand(Float32, dim, m))
    optqueries = StrideMatrixDatabase(rand(Float32, dim, 10))
    E = ExhaustiveSearch(; dist, db)
    bruteforcesearchtime = @elapsed Igold, Dgold = searchbatch(E, queries, k)
    minleaf = 100
    numqueries = 16  # numqueries for optimize!
    verbose = false

    for e in [
            (sortsat=ProximalSortSat(), exact=0.9999, arecall=0.8, brecall=0.3, crecall=0.5),
            (sortsat=RandomSortSat(), exact=0.9999, arecall=0.8, brecall=0.6, crecall=0.6),
            (sortsat=DistalSortSat(), exact=0.9999, arecall=0.8, brecall=0.6, crecall=0.6)
            ]
        sortsat = e.sortsat
        @info "============================= $sortsat minleaf=$minleaf ================"
        @info "Sat"
        sat = Sat(db; dist)
        @time index!(sat; sortsat, minleaf)
        # @show sat
        # checking that the database size and the number of inserted elements is consistent
        @test n == 1 + sum(length(C) for C in sat.children if C !== nothing)
        # checking that cov is consistent with children
        #@test all((C === nothing ? sat.cov[i] < 0 : sat.cov[i] > 0) for (i, C) in enumerate(sat.children))
        @test all(sat.cov[i] >= 0  for (i, C) in enumerate(sat.children))
        Isat, Dsat = searchbatch(sat, queries, k)
        searchtime = @elapsed Isat, Dsat = searchbatch(sat, queries, k)
        recall = macrorecall(Igold, Isat)
        @test recall >= e.exact
        continue
        
        asat = PruningSat(sat)
        alist = optimize!(asat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        Ia, _ = searchbatch(asat, queries, k)
        asearchtime = @elapsed Ia, _ = searchbatch(asat, queries, k)
        arecall = macrorecall(Igold, Ia)
        @test arecall > e.arecall ## it should be close to 0.9, but anyway errors will show very low recall

        bsat = BeamSearchSat(sat)
        blist = optimize!(bsat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        #@show bsat.bs
        Ib, _ = searchbatch(bsat, queries, k)
        bsearchtime = @elapsed Ib, _ = searchbatch(bsat, queries, k)
        brecall = macrorecall(Igold, Ib)
        @test brecall > e.brecall
        @info "multisat"

        @time csat = BeamSearchMultiSat([index!(Sat(db; dist, root=rand(1:n)); sortsat, minleaf) for _ in 1:4])
        # @show [s.root for s in csat.sat]
        clist = optimize!(csat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        #@show csat.bs
        Ic, _ = searchbatch(csat, queries, k)
        csearchtime = @elapsed Ic, _ = searchbatch(csat, queries, k)
        crecall = macrorecall(Igold, Ic)
        @show crecall > e.crecall

        @info "SearchGraph"
        @time G = index!(SearchGraph(; db, dist, verbose))
        # @show [s.root for s in csat.sat]
        glist = optimize!(G, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        #@show csat.bs
        Ig, _ = searchbatch(G, queries, k)
        gsearchtime = @elapsed Ig, _ = searchbatch(G, queries, k)
        grecall = macrorecall(Igold, Ig)
        @test grecall > 0.6

        @info "------------ brute force: $bruteforcesearchtime -- sort sat: $(typeof(sortsat)) "
        @info " exact sat:" recall searchtime (bruteforcesearchtime / searchtime)
        @info " sat with probabilistic spell:" arecall asearchtime (bruteforcesearchtime / asearchtime) first(alist)
        @info " sat with beam search:" brecall bsearchtime (bruteforcesearchtime / bsearchtime) first(blist)
        @info " multi sat with beam search:" crecall csearchtime (bruteforcesearchtime / csearchtime) first(clist)
        @info " search graph:" grecall gsearchtime (bruteforcesearchtime / gsearchtime) first(glist)
    end
end

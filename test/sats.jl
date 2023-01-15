# this file is part of SpatialAccessTrees.jl

using SpatialAccessTrees, SimilaritySearch
using Test

@testset "SpatialAccessTrees" begin
    n, m, k, dim, minleaf, dist = 10^4, 10^2, 10, 4, 4, L2Distance()
    db = StrideMatrixDatabase(rand(Float32, dim, n))
    warmqueries = StrideMatrixDatabase(rand(Float32, dim, 2))
    queries = StrideMatrixDatabase(rand(Float32, dim, m))
    optqueries = StrideMatrixDatabase(rand(Float32, dim, 10))
    E = ExhaustiveSearch(; dist, db)
    Igold, _ = searchbatch(E, warmqueries, k)
    bruteforcesearchtime = @elapsed Igold, _ = searchbatch(E, queries, k)
    numqueries = 16  # numqueries for optimize!
    verbose = false

    for e in [
            (ipart=SatInitialPartition(), sortsat=ProximalSortSat(), exact=0.9999, arecall=0.75, brecall=0.3, crecall=0.5),
            (ipart=SatInitialPartition(), sortsat=DistalSortSat(), exact=0.9999, arecall=0.75, brecall=0.6, crecall=0.6),
            # it seems that we destroy some properties used by BeamSearchSat using RandomInitialPartition
            (ipart=RandomInitialPartition(nparts=64, shuffle=false), sortsat=RandomSortSat(), exact=0.9999, arecall=0.75, brecall=0.1, crecall=0.6),
            (ipart=RandomInitialPartition(nparts=64, shuffle=true), sortsat=RandomSortSat(), exact=0.9999, arecall=0.75, brecall=0.1, crecall=0.6),
            (ipart=SatInitialPartition(), sortsat=RandomSortSat(), exact=0.9999, arecall=0.75, brecall=0.6, crecall=0.6),
        ]
        sortsat, ipart = e.sortsat, e.ipart
        @info "============================= $ipart -- $sortsat -- minleaf=$minleaf ================"
        @info "Sat"
        sat = Sat(db; dist)
        buildtime = @elapsed index!(sat, ipart; sortsat, minleaf)
        @show sum(c !== nothing for c in sat.children)
        # @show sat
        # checking that the database size and the number of inserted elements is consistent
        @test n == 1 + sum(length(C) for C in sat.children if C !== nothing)
        # checking that cov is consistent with children
        #@test all((C === nothing ? sat.cov[i] < 0 : sat.cov[i] > 0) for (i, C) in enumerate(sat.children))
        @test all(sat.cov[i] >= 0  for (i, C) in enumerate(sat.children))
        Isat, _ = searchbatch(sat, warmqueries, k)
        searchtime = @elapsed Isat, _ = searchbatch(sat, queries, k)
        recall = macrorecall(Igold, Isat)
        @test recall >= e.exact

        psat = permutesat(sat)
        @test isperm(psat.π)
        Ip, _ = searchbatch(psat, warmqueries, k)
        psearchtime = @elapsed Ip, _ = searchbatch(psat, queries, k)
        precall = macrorecall(Igold, Ip)
        @test precall >= e.exact
        
        asat = PruningSat(sat)
        alist = optimize!(asat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        ## asat = PruningSat(psat.index)
        ## alist = optimize!(asat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        ## asat = PermutedSearchIndex(index=asat, π=psat.π, π′=psat.π′)
        Ia, _ = searchbatch(asat, warmqueries, k)
        asearchtime = @elapsed Ia, _ = searchbatch(asat, queries, k)
        arecall = macrorecall(Igold, Ia)
        @test arecall > e.arecall ## it should be close to 0.9, but anyway errors will show very low recall

        pasat = PruningSat(psat.index)
        palist = optimize!(pasat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        pasat = PermutedSearchIndex(index=pasat, π=psat.π, π′=psat.π′)
        Ipa, _ = searchbatch(pasat, warmqueries, k)
        pasearchtime = @elapsed Ipa, _ = searchbatch(pasat, queries, k)
        parecall = macrorecall(Igold, Ipa)
        @test parecall > e.arecall ## it should be close to 0.9, but anyway errors will show very low recall

        bsat = BeamSearchSat(sat)
        blist = optimize!(bsat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        # @show bsat.bs
        Ib, _ = searchbatch(bsat, warmqueries, k)
        bsearchtime = @elapsed Ib, _ = searchbatch(bsat, queries, k)
        brecall = macrorecall(Igold, Ib)
        @test brecall > e.brecall
        
        #=
        @info "multisat"

        @time csat = BeamSearchMultiSat([index!(Sat(db; dist, root=rand(1:n)); sortsat, minleaf) for _ in 1:4])
        # @show [s.root for s in csat.sat]
        clist = optimize!(csat, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        #@show csat.bs
        Ic, _ = searchbatch(csat, queries, k)
        csearchtime = @elapsed Ic, _ = searchbatch(csat, queries, k)
        crecall = macrorecall(Igold, Ic)
        @show crecall > e.crecall
        =#
        @info "SearchGraph"
        @time G = index!(SearchGraph(; db, dist, verbose))
        # @show [s.root for s in csat.sat]
        glist = optimize!(G, MinRecall(0.9); numqueries, verbose, ksearch=k, queries=optqueries)
        #@show csat.bs
        Ig, _ = searchbatch(G, queries, k)
        gsearchtime = @elapsed Ig, _ = searchbatch(G, queries, k)
        grecall = macrorecall(Igold, Ig)
        @test grecall > 0.6

        @info "------------ brute force: $bruteforcesearchtime -- sort sat: $(typeof(sortsat)) -- brute-force searchtime: $bruteforcesearchtime, buildtime: $buildtime"
        @info " exact sat:" recall searchtime (bruteforcesearchtime / searchtime)
        @info " permuted exact sat:" precall psearchtime (bruteforcesearchtime / psearchtime)
        @info " sat with probabilistic spell:" arecall asearchtime (bruteforcesearchtime / asearchtime) first(alist)
        @info " perm sat with probabilistic spell:" parecall pasearchtime (bruteforcesearchtime / pasearchtime) first(palist)
        @info " sat with beam search:" brecall bsearchtime (bruteforcesearchtime / bsearchtime) first(blist)
        #@info " multi sat with beam search:" crecall csearchtime (bruteforcesearchtime / csearchtime) first(clist)
        @info " search graph:" grecall gsearchtime (bruteforcesearchtime / gsearchtime) first(glist)
    end
end

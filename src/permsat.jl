# this file is part of SpatialAccessTrees.jl

export permutesat

function naivesatpermutation!(π, sat::Sat)
    p = 1
    π[p] = sat.root
    
    for C in sat.children
        C === nothing && continue
        for i in C
            p += 1
            π[p] = i
        end
    end

    π
end

function satpermutation!(π, sat::Sat)
    p = 1
    π[p] = sat.root
    cand = [sat.root]
    sizehint!(cand, ceil(Int, sqrt(length(π))))

    while length(cand) > 0
        i = popfirst!(cand)
        for C in sat.children[i]
            for c in C
                p += 1
                π[p] = c
                if sat.children[c] !== nothing
                    push!(cand, c)
                end
            end
        end
    end

    π
end

satpermutation(sat::Sat) = satpermutation!(Vector{UInt32}(undef, length(sat)), sat)

"""
    permutesat(sat::Sat, π=satpermutation(sat), π′=invperm(π))

Permute sat to optimize cache accesses patterns; the database is also copied and permuted.
The permuted index is stored in a `PermutedSearchIndex` struct to allow plug and play index interchange.

# Arguments:
- `sat`: Input `Sat` index.
- `π`: Permutation.
- `π′`: Inverse permutation.
"""
function permutesat(sat::Sat, π=satpermutation(sat), π′=invperm(π))
    db = MatrixDatabase(SubDatabase(database(sat), π))
    children = similar(sat.children)
    fill!(children, nothing)

    for ii in eachindex(children)
        C = sat.children[ii]
        C === nothing && continue
        CC = copy(C)
        children[π′[ii]] = CC

        for (i, c) in enumerate(C)
            CC[i] = π′[c]
        end
    end

    cov = similar(sat.cov)
    for i in eachindex(cov)
        cov[π′[i]] = sat.cov[i]
    end

    s = Sat(distance(sat), db, π[sat.root], children, cov)
    PermutedSearchIndex(s, π, π′)
end

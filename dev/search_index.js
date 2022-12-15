var documenterSearchIndex = {"docs":
[{"location":"","page":"API","title":"API","text":"\nCurrentModule = SpatialAccessTrees\nDocTestSetup = quote\n    using SpatialAccessTrees\nend","category":"page"},{"location":"#Indexes","page":"API","title":"Indexes","text":"","category":"section"},{"location":"","page":"API","title":"API","text":"Sat\nindex!","category":"page"},{"location":"#SpatialAccessTrees.Sat","page":"API","title":"SpatialAccessTrees.Sat","text":"struct Sat{DT<:SemiMetric,DBT<:AbstractDatabase} <: AbstractSearchIndex\n    dist::DT\n    db::DBT\n    root::UInt32\n    parents::Vector{UInt32}\n    children::Vector{Union{Nothing,Vector{UInt32}}}\n    cov::Vector{Float32} # leafs: ``- d(parent, leaf)``, internal: ``\\max \\{d(parent, u) | u \\in children(parent)\\}``\nend\n\nSpatial Access Tree data structure. Please see Sat function for high level constructors\n\n\n\n\n\n","category":"type"},{"location":"#SimilaritySearch.index!","page":"API","title":"SimilaritySearch.index!","text":"index!(\n    sat::Sat;\n    sortsat::AbstractSortSat=ProximalSortSat(),\n    minleaf::Int=log2(ceil(database(sat)))\n)\n\nPerforms the indexing of the referenced dataset in the spatial access tree.\n\nArguments\n\nsat: The metric data structure\n\nKeyword arguments\n\nsortsat: The strategy to create the spatial access tree, it heavily depends on the order of elements while it is build. It accepts:\nRandomSortSat(): children are randomized.\nProximalSortSat(): classical approach, near elements are put first.\nDistalSortSat(): recent approach, distant elements are put first.\nminleaf: Minimum number of children to perform a spatial access separation (half space partitioning)\n\n\n\n\n\n","category":"function"},{"location":"#Searching","page":"API","title":"Searching","text":"","category":"section"},{"location":"","page":"API","title":"API","text":"search\n","category":"page"},{"location":"","page":"API","title":"API","text":"@autodocs Modules = [SpatialAccessTrees] Public = false Order = [:function, :type] ```","category":"page"}]
}
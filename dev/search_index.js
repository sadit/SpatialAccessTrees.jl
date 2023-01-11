var documenterSearchIndex = {"docs":
[{"location":"","page":"API","title":"API","text":"\nCurrentModule = SpatialAccessTrees\nDocTestSetup = quote\n    using SpatialAccessTrees\nend","category":"page"},{"location":"#Spatial-Access-Trees-(SAT)","page":"API","title":"Spatial Access Trees (SAT)","text":"","category":"section"},{"location":"","page":"API","title":"API","text":"Sat\nindex!","category":"page"},{"location":"#SpatialAccessTrees.Sat","page":"API","title":"SpatialAccessTrees.Sat","text":"struct Sat{DT<:SemiMetric,DBT<:AbstractDatabase} <: AbstractSearchIndex\n    dist::DT\n    db::DBT\n    root::UInt32\n    # parents::Vector{UInt32}\n    children::Vector{Union{Nothing,Vector{UInt32}}}\n    cov::Vector{Float32} # leafs: ``- d(parent, leaf)``, internal: ``\\max \\{d(parent, u) | u \\in children(parent)\\}``\nend\n\nSpatial Access Tree data structure. Please see Sat function for high level constructors\n\n\n\n\n\n","category":"type"},{"location":"#SimilaritySearch.index!","page":"API","title":"SimilaritySearch.index!","text":"index!(sat::Sat; <kwargs...>)\nindex!(sat::Sat, ipart::SatInitialPartition; <kwargs...>)\nindex!(sat::Sat, ipart::RandomInitialPartition; <kwargs...>)\n\nPerforms the indexing of the referenced dataset in the tree. It supports limited forms of multithreading, induced by initial partitioning schemes.\n\nArguments\n\nsat: The metric data structure.\nipart: initial partitioning scheme for the tree. It supports the following kinds of objects:\nSatInitialPartition(): Traditional construction, default value. Each part is a SAT partition and will be processed by a thread on multithreading setups.\nRandomInitialPartition(nparts=Threads.nthreads(), shuffle=false):   construction that divides the dataset (randomly if suffle=true) in nparts disjoint parts. The resulting structure violates the SAT partitioning in a whole   and creates a kind of SAT forest that are fine SAT partitions. Useful to limit the height of the tree and for multiprocessing purporses, i.e., each part will be processed by a thread.\n\nKeyword arguments\n\nsortsat: The strategy to create the spatial access tree, it heavily depends on the order of elements while it is build. It accepts:\nRandomSortSat(): children are randomized (default value)\nProximalSortSat(): classical approach, near elements are put first.\nDistalSortSat(): recent approach, distant elements are put first.\nminleaf: Minimum number of children to perform a spatial access separation (half space partitioning)\n\n\n\n\n\n","category":"function"},{"location":"#SAT-orders","page":"API","title":"SAT orders","text":"","category":"section"},{"location":"","page":"API","title":"API","text":"ProximalSortSat\nDistalSortSat\nRandomSortSat","category":"page"},{"location":"#Approximate-similarity-search","page":"API","title":"Approximate similarity search","text":"","category":"section"},{"location":"","page":"API","title":"API","text":"PruningSat\nBeamSearchSat\nBeamSearchMultiSat\n\noptimize!","category":"page"},{"location":"#SpatialAccessTrees.PruningSat","page":"API","title":"SpatialAccessTrees.PruningSat","text":"PruningSat(sat::Sat; factor=0.9)\n\nCreates an approximate similarity sarch index based on sat and aggressive pruning; adapted from the paper A probabilistic spell for the curse of dimensionality (Chavez and Navarro, 2001). It supports auto-tuning via  optimize!.\n\n\n\n\n\n","category":"type"},{"location":"#SpatialAccessTrees.BeamSearchSat","page":"API","title":"SpatialAccessTrees.BeamSearchSat","text":"BeamSearchSat(sat::Sat; bsize=8, Δ=1f0)\n\nCreates an approximate similarity sarch index based on sat and aggressive pruning; adapted from the paper A probabilistic spell for the curse of dimensionality (Chavez and Navarro, 2001). It supports auto-tuning via  optimize!.\n\n\n\n\n\n","category":"type"},{"location":"#SpatialAccessTrees.BeamSearchMultiSat","page":"API","title":"SpatialAccessTrees.BeamSearchMultiSat","text":"BeamSearchMultiSat(sat::Vector{SatType}; bsize=8, Δ=1f0) where {SatType<:Sat}\n\nApplies beam search on a multi SAT array. See optimize! to tune the index for some objective.\n\n\n\n\n\n","category":"type"},{"location":"","page":"API","title":"API","text":"Modules = [SpatialAccessTrees]\nPublic = false\nOrder = [:function, :type]","category":"page"}]
}

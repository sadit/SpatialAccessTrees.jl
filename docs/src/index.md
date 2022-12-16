```@meta

CurrentModule = SpatialAccessTrees
DocTestSetup = quote
    using SpatialAccessTrees
end
```

# Spatial Access Trees (SAT)
```@docs
Sat
index!
```

## SAT orders

```@docs
ProximalSortSat
DistalSortSat
RandomSortSat
```

# Approximate similarity search
```@docs
PruningSat
BeamSearchSat
BeamSearchMultiSat

optimize!
```

```@autodocs
Modules = [SpatialAccessTrees]
Public = false
Order = [:function, :type]
```

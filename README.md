# SpatialAccessTrees

Spatial access trees are a family of metric trees having excellent performance on low and medium dimensional datasets. 

This Julia package implements some of metric indexes based on the spatial access tree described in:

```
Navarro, G. (2002). Searching in metric spaces by spatial approximation. The VLDB Journal, 11(1), 28-46.
```

and

```
Chávez, E., Luduena, V., Reyes, N., & Roggero, P. (2016). Faster proximity searching with the distal SAT. Information Systems, 59, 15-47.
```

The package supports for approximate similarity search using aggressive pruning

```
Chávez, E., & Navarro, G. (2001, January). A probabilistic spell for the curse of dimensionality. In Workshop on Algorithm Engineering and Experimentation (pp. 147-160). Springer, Berlin, Heidelberg.
```


I will also add some new structures related to them.


This package was designed to work with [`SimilaritySearch.jl`](https://github.com/sadit/SimilaritySearch.jl).

using Documenter, SpatialAccessTrees

makedocs(;
    modules=[SpatialAccessTrees],
    authors="Eric S. Tellez",
    repo="https://github.com/sadit/SpatialAccessTrees.jl/blob/{commit}{path}#L{line}",
    sitename="SpatialAccessTrees.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://sadit.github.io/SpatialAccessTrees.jl",
        assets=String[],
    ),
    pages=[
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/sadit/SpatialAccessTrees.jl",
    devbranch=nothing,
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)

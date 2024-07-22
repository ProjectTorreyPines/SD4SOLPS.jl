using Documenter
using SD4SOLPS

makedocs(;
    modules=[SD4SOLPS],
    format=Documenter.HTML(),
    sitename="SD4SOLPS",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/SD4SOLPS.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"],
)

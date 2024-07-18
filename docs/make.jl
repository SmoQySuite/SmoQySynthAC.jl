using SmoQySynthAC
using Documenter
using Literate

DocMeta.setdocmeta!(SmoQySynthAC, :DocTestSetup, :(using SmoQySynthAC); recursive=true)

# package directory
pkg_dir = pkgdir(SmoQySynthAC)
# docs directory
docs_dir = joinpath(pkg_dir, "docs")
# docs/src directory
docs_src_dir = joinpath(docs_dir, "src")

# generate usage markdown file
usage_jl = joinpath(docs_dir, "usage.jl")
Literate.markdown(usage_jl, docs_src_dir; credit = false)

makedocs(;
    clean = true,
    modules=[SmoQySynthAC],
    authors="Benjamin Cohen-Stead <benwcs@gmail.com>, Steven Johnston <sjohn145@utk.edu>",
    sitename="SmoQySynthAC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SmoQySuite.github.io/SmoQySynthAC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Usage" => "usage.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/SmoQySuite/SmoQySynthAC.jl",
    devbranch="main",
)

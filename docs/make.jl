using SmoQySynthAC
using Documenter

DocMeta.setdocmeta!(SmoQySynthAC, :DocTestSetup, :(using SmoQySynthAC); recursive=true)

makedocs(;
    modules=[SmoQySynthAC],
    authors="Benjamin Cohen-Stead <benwcs@gmail.com>, Steven Johnston <sjohn145@utk.edu>",
    sitename="SmoQySynthAC.jl",
    format=Documenter.HTML(;
        canonical="https://SmoQySuite.github.io/SmoQySynthAC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SmoQySuite/SmoQySynthAC.jl",
    devbranch="main",
)

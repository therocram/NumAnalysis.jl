# include("../src/NumAnalysis.jl")
# using Documenter, .NumAnalysis

using Documenter, NumAnalysis

DocMeta.setdocmeta!(NumAnalysis, :DocTestSetup, :(using NumAnalysis))

makedocs(
    modules = [NumAnalysis],

    pages = [
        "Home" => "index.md",
        "Numerical Integration" => "numintegration.md",
        "Initial Value Problems" => "initvalode.md",
        "Root Finder" => "rootfinder.md"
    ],

    sitename="NumAnalysis.jl"
)
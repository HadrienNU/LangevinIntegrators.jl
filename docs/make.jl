using LangevinIntegrators
using LangevinIntegratorsPlumedExt
using Documenter

DocMeta.setdocmeta!(LangevinIntegrators, :DocTestSetup, :(using LangevinIntegrators); recursive=true)

makedocs(;
    modules=[LangevinIntegrators],
    authors="Hadrien <hadrien.vroylandt@sorbonne-universite.fr> and contributors",
    repo="https://github.com/HadrienNU/LangevinIntegrators.jl/blob/{commit}{path}#{line}",
    sitename="LangevinIntegrators.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HadrienNU.github.io/LangevinIntegrators.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Integrators" => "integrators.md",
        "Forces and fixes" => "forces.md",
        "Plumed" => "plumed.md",
    ],
)

deploydocs(;
    repo="github.com/HadrienNU/LangevinIntegrators.jl",
    devbranch="integrators",
)

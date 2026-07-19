using AbstractQAtlas
using Documenter
using Downloads

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)
favicon_path = joinpath(assets_dir, "favicon.ico")
logo_path = joinpath(assets_dir, "logo.png")

Downloads.download("https://github.com/sotashimozono.png", favicon_path)
Downloads.download("https://github.com/sotashimozono.png", logo_path)

makedocs(;
    sitename="AbstractQAtlas.jl",
    format=Documenter.HTML(;
        canonical="https://codes.sota-shimozono.com/AbstractQAtlas.jl/stable/",
        prettyurls=get(ENV, "CI", "false") == "true",
        mathengine=MathJax3(
            Dict(
                :tex => Dict(
                    :inlineMath => [["\$", "\$"], ["\\(", "\\)"]],
                    :tags => "ams",
                    :packages => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
        assets=["assets/favicon.ico", "assets/custom.css"],
        # The single @autodocs API page grows as relations accumulate; lift
        # the HTML size gate for now.  TODO(refactor): split the API
        # reference into per-domain pages once it settles.
        size_threshold=1_000_000,
        size_threshold_warn=1_000_000,
    ),
    modules=[AbstractQAtlas],
    pages=["Home" => "index.md"],
)

deploydocs(;
    versions=["stable", "dev"],
    repo="github.com/QAtlasHub/AbstractQAtlas.jl.git",
    devbranch="main",
    push_preview=true,
)

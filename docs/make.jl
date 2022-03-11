using Documenter
using TFT

makedocs(
    modules     = [TFT],
    format      = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename    = "TFT.jl",
    authors     = "Tom Van Acker",
    pages       = [ "Home"              => "index.md",
                  ]
)

deploydocs(
     repo = "github.com/timmyfaraday/TFT.jl.git"
)
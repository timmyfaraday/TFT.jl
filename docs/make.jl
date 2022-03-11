using Documenter
using TFT

makedocs(
    sitename = "TFT",
    format = Documenter.HTML(),
    modules = [TFT]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#



# local
using Pkg
Pkg.activate(".")

using Documenter, Revise, VortexDistributions
makedocs(sitename="VortexDistributions Documentation",
format = Documenter.HTML(prettyurls = false))

# CI
# using Documenter, VortexDistributions
# makedocs(sitename="VortexDistributions Documentation",
#     format = Documenter.HTML(
#     prettyurls = get(ENV, "CI", nothing) == "true")
#     )

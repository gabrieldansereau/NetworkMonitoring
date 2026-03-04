# using DrWatson
# @quickactivate :NetworkMonitoring

# Note: the NetworkMonitoring module cannot be loaded at the same time as
# AlgebraOfGraphics due to a weird name conflicts. As a temporary workaround, we
# can load only the packages needed for plots and summary analyses, and
# includet() to track certain files from the module with Revise.

using AlgebraOfGraphics
using CairoMakie
using Combinatorics
using CSV
using DataFramesMeta
using DrWatson
using ProgressMeter
using Random
using Revise
using Statistics
import SpeciesDistributionToolkit as SDT

update_theme!(; CairoMakie=(; px_per_unit=2.0))
CairoMakie.activate!(; type="svg")

includet(srcdir("colours.jl"))
includet(srcdir("efficiency.jl"))
includet(srcdir("focal.jl"))
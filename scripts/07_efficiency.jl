# using DrWatson
# @quickactivate :NetworkMonitoring

using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFramesMeta
using DrWatson
using Statistics
import SpeciesDistributionToolkit as SDT

update_theme!(; CairoMakie=(; px_per_unit=2.0))

## Load simulations

# Load summarized simulation results
sims_samplers = CSV.read(datadir("sims_samplers.csv"), DataFrame)
sims_optimized = CSV.read(datadir("sims_optimized.csv"), DataFrame)

# Select UncertaintySampling only as example
sims_set = @rsubset(sims_samplers, :sampler == "UncertaintySampling")

# Visualize
begin
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites", ylabel="Proportion")
    for r in unique(sims_set.sim)
        _sim = @rsubset(sims_set, :sim == r)
        scatter!(ax, _sim.nbon, _sim.med; color=Makie.wong_colors()[2])
    end
    fig
end

## Efficiency

# Saturation equation
saturation(a) = (x) -> x ./ (a .+ x)

# Efficiency grid search
function efficiency(x, y; A=LinRange(-12.0, 12.0, 500))
    err = zeros(length(A))
    for i in eachindex(A)
        f = saturation(exp(A[i]))
        err[i] = sqrt(sum((y .- f(x)) .^ 2.0))
    end
    return exp(A[last(findmin(err))])
end

# Visualize
begin
    f = Figure()
    ax = Axis(f[1, 1])
    for iter in unique(sims_set.sim)
        X₁ = @rsubset(sims_set, :sim == iter)
        x = X₁.nbon
        y = X₁.med
        eff = efficiency(x, y)
        scatter!(ax, x, y; color=:lightgrey)
        lines!(ax, x, saturation(eff)(x); color=:black, linestyle=:dash)
    end
    f
end

# Example
sim1 = @rsubset(sims_set, :sim == 1)
x = sim1.nbon
y = sim1.med
eff = efficiency(x, y)
saturation(eff)(x)

## Occupancy

# Load layer for exploration
_l = SDT.SDMLayer(datadir("focal_array", "layer_sp_range-01.tiff"))
heatmap(_l)

# Get layer occupancy
occupancy(_l) = length(findall(isone, _l)) / length(_l)
occupancy(_l)

# Occupancy across multiple layers
occup = zeros(10)
for i in 1:10
    idp = lpad(i, 2, "0")
    l = SDT.SDMLayer(datadir("focal_array", "layer_sp_range-$idp.tiff"))
    occup[i] = occupancy(l)
end
occup

# Corresponding sampling efficiency
effs = zeros(10)
for i in 1:10
    X₁ = @rsubset(sims_set, :sim == i)
    x = X₁.nbon
    y = X₁.med
    effs[i] = efficiency(x, y)
end
effs

# Visualize
scatter(occup, effs)

# AOG
mapping(occup, effs) * visual(Scatter) |> draw
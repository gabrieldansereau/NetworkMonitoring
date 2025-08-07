# using DrWatson
# @quickactivate :NetworkMonitoring

using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFramesMeta
using DrWatson
using Random
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
mapping(occup, effs) * AlgebraOfGraphics.density() |> draw
mapping(occup, effs) * AlgebraOfGraphics.density() * visual(Contour) |> draw

## Scale up comparison

# Calculate efficiency & assign occupancy
effs_samplers = @chain sims_samplers begin
    @groupby(:sim, :sampler)
    @combine(:eff = efficiency(:nbon, :med))
    @transform(:occ = occup[:sim], :set = "Samplers")
end
effs_optimized = @chain sims_optimized begin
    @groupby(:sim, :sampler)
    @combine(:eff = efficiency(:nbon, :med))
    @transform(:occ = occup[:sim], :set = "Layers")
end

# Define color sets
cols = [
    # Interaction types
    "possible" => Makie.wong_colors()[2],
    "realized" => Makie.wong_colors()[3],
    "detected" => Makie.wong_colors()[4],
    # Samplers
    "UncertaintySampling" => Makie.wong_colors()[2],
    "WeightedBalancedAcceptance" => Makie.wong_colors()[3],
    "SimpleRandom" => Makie.wong_colors()[1],
    # Layers
    "Focal species range" => Makie.wong_colors()[2],
    "Species richness" => Makie.wong_colors()[4],
    "Realized interactions" => Makie.wong_colors()[5],
]

# Scatter & smooth
efflog = :eff => log => "log(eff)"
layout = mapping(:occ, efflog; color=:sampler) * (visual(Scatter) + smooth())
fig =
    data(effs_samplers) * layout |>
    x -> draw(x, scales(; Color=(; palette=cols)); legend=(; position=:bottom))
fig =
    data(effs_optimized) * layout |>
    x -> draw(x, scales(; Color=(; palette=cols)); legend=(; position=:bottom))

# Heatmaps
layout =
    mapping(:occ, efflog; row=:sampler) * (AlgebraOfGraphics.density() + visual(Scatter))
f1 = data(effs_samplers) * layout
f2 = data(effs_optimized) * layout
draw(f1)
draw(f2)
begin
    f = Figure()
    draw!(f[1, 1], f1)
    draw!(f[1, 2], f2)
    f
end

# Violin
begin
    Random.seed!(42) # for jitter
    layer =
        mapping(:sampler, efflog; color=:sampler) *
        visual(RainClouds; markersize=10, jitter_width=0.1, plot_boxplots=false)
    f1 = data(effs_samplers) * layer
    f2 = data(effs_optimized) * layer
    f = Figure()
    draw!(f[1, 1], f1, scales(; Color=(; palette=cols)))
    draw!(f[2, 1], f2, scales(; Color=(; palette=cols)))
    f
end

## Within-simulation variation

# Load one simulation example
monitored_samplers = CSV.read(datadir("monitored_samplers.csv"), DataFrame)
@rsubset!(monitored_samplers, :sampler == "UncertaintySampling")

# Visualize
let b = monitored_samplers
    eff = efficiency(b.nbon, b.med)
    band(b.nbon, b.low, b.upp; alpha=0.4, color=Makie.wong_colors()[2])
    scatter!(b.nbon, b.med; color=Makie.wong_colors()[2])
    lines!(b.nbon, saturation(eff)(b.nbon); color=:black, linestyle=:dash)
    current_figure()
end

# Load pre-summary results
sim1_samplers = @chain begin
    CSV.read(datadir("focal_array", "monitored_samplers-01.csv"), DataFrame)
    @rsubset(:sampler == "UncertaintySampling")
    @rtransform(:monitored = :monitored / :deg)
end

# Visualize
fig = let bs = sim1_samplers, b = monitored_samplers
    f = Figure()
    ax = Axis(f[1, 1]; xlabel="Number of sites", ylabel="Monitored proportion")
    for i in 1:100
        bi = @rsubset(bs, :rep == i)
        # lines!(bi.nbon, bi.monitored; color=:grey, linewidth=0.5)
        eff = efficiency(bi.nbon, bi.monitored; A=LinRange(-12.0, 12.0, 5000))
        lines!(
            bi.nbon,
            saturation(eff)(b.nbon);
            color=:grey,
            linestyle=:solid,
            label="individual saturation curves",
        )
    end
    eff = efficiency(b.nbon, b.med)
    band!(
        b.nbon,
        b.low,
        b.upp;
        alpha=0.4,
        color=Makie.wong_colors()[2],
        label="90% percentile",
    )
    scatter!(b.nbon, b.med; color=Makie.wong_colors()[2], label="median points")
    lines!(
        b.nbon,
        saturation(eff)(b.nbon);
        color=:black,
        linestyle=:dash,
        label="saturation from median",
    )
    axislegend(; position=:rb, unique=true)
    current_figure()
end
save(plotsdir("saturation_within_example.png"), fig)

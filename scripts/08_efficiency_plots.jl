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

# Load data
sims_samplers = CSV.read(datadir("sims_efficiency_samplers.csv"), DataFrame)
sims_optimized = CSV.read(datadir("sims_efficiency_optimized.csv"), DataFrame)
sims_species = CSV.read(datadir("sims_efficiency_species.csv"), DataFrame)
effs_samplers = CSV.read(datadir("efficiency_samplers.csv"), DataFrame)
effs_optimized = CSV.read(datadir("efficiency_optimized.csv"), DataFrame)
effs_species = CSV.read(datadir("efficiency_species.csv"), DataFrame)

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

## Efficiency only

# Violin
begin
    Random.seed!(42) # for jitter
    efflog = :eff => log => "log(efficiency)"
    layer =
        mapping(:sampler, efflog; color=:sampler) *
        visual(RainClouds; markersize=5, jitter_width=0.1, plot_boxplots=false)
    f1 = data(effs_samplers) * layer
    f2 = data(effs_optimized) * layer
    f = Figure(; size=(700, 700))
    draw!(f[1, 1], f1, scales(; Color=(; palette=cols)))
    draw!(f[2, 1], f2, scales(; Color=(; palette=cols)))
    f
end
save(plotsdir("efficiency_distribution.png"), f)

# Efficiency given species degree
begin
    Random.seed!(42)
    ranks = :rank => ranknames => "Within-simulation Degree Percentile Rank"
    f = Figure(; size=(650, 650))
    fg1 = draw!(
        f[1, 1:2],
        data(effs_species) *
        mapping(ranks, efflog; color=ranks) *
        visual(RainClouds; markersize=5, jitter_width=0.1, plot_boxplots=false),
        scales(; Color=(; palette=from_continuous(cgrad(:viridis; rev=true)))),
    )
    fg2 = draw!(
        f[2, 1],
        data(effs_species) *
        mapping(:deg => "Degree", efflog; color=:rank => ranknames => "Percentile rank") *
        visual(Scatter),
        scales(; Color=(; palette=from_continuous(cgrad(:viridis; rev=true)))),
    )
    legend!(f[2, 2], fg2)
    f
end
save(plotsdir("efficiency_distribution_species.png"), f)

## Efficiency & Occupancy

# Scatter & smooth
sortedsamplers = ["UncertaintySampling", "WeightedBalancedAcceptance", "SimpleRandom"]
sortedlayers = ["Focal species range", "Species richness", "Realized interactions"]
sortedlayout = [sortedsamplers..., sortedlayers...]
begin
    occ = :occ => "occupancy"
    layout = mapping(occ, efflog; color=:sampler) * (visual(Scatter) + linear())
    scale = scales(; Color=(; palette=cols))
    legend = (; position=:bottom)
    f1 = data(effs_samplers) * layout
    f2 = data(effs_optimized) * layout * mapping(; color=:sampler => "optimization layer")
end
draw(f1 * mapping(; col=:sampler => sorter(sortedsamplers)), scale; legend=legend)
draw(f2 * mapping(; col=:sampler => sorter(sortedlayers)), scale; legend=legend)
fig = draw(
    (f1 + f2) * mapping(; layout=:sampler => sorter(sortedlayout)),
    scales(; Color=(; palette=cols, legend=false));
    figure=(; size=(800, 450)),
)
save(plotsdir("efficiency_occupancy.png"), fig)

# Same with independent subfigures
fig = let
    f = Figure(; size=(800, 450))
    fg1 = draw!(f[1, 1], f1 * mapping(; col=:sampler => sorter(sortedsamplers)), scale)
    fg2 = draw!(f[2, 1], f2 * mapping(; col=:sampler => sorter(sortedlayers)), scale)
    linkyaxes!(fg1..., fg2...)
    f
end

# Group columns in single panel
fig = let
    f = Figure(; size=(800, 450))
    fg1 = draw!(f[1, 1], f1, scale)
    legend!(f[2, 1], fg1; position=:bottom, tellheight=true, tellwidth=false)
    fg2 = draw!(f[1, 2], f2, scale)
    legend!(f[2, 2], fg2; position=:bottom, tellheight=true, tellwidth=false)
    linkyaxes!(fg1..., fg2...)
    f
end
save(plotsdir("_xtras/", "efficiency_occupancy_scatter.png"), fig)

# Species degree & efficiency-occupancy
fig = draw(data(effs_species) * layout * mapping(; color=:deg); legend=legend)
save(plotsdir("_xtras/", "efficiency_occupancy_species_degree.png"), fig)

# Species rank instead of degree
ranknames = renamer(1 => "1.0", 2 => "0.66", 3 => "0.33", 4 => "0.07")
fig = draw(
    data(effs_species) *
    layout *
    mapping(; color=:rank => ranknames => "Within-simulation Degree Percentile Rank"),
    scales(; Color=(; palette=from_continuous(cgrad(:viridis; rev=true))));
    legend=legend,
)
save(plotsdir("_xtras", "efficiency_occupancy_species_rank.png"), fig)

# Same with facets
fig = draw(
    data(effs_species) * layout * mapping(; color=ranks, layout=ranks),
    scales(; Color=(; palette=from_continuous(cgrad(:viridis; rev=true))));
    figure=(; size=(500, 450)),
    legend=legend,
)
save(plotsdir("efficiency_occupancy_degree_ranked_facets.png"), fig)

## Within-simulation comparison

# NDI
ndi(x, y) = (x - y) / (x + y)

# Separate results per simulation
effs_combined = vcat(effs_samplers, effs_optimized)
function comparewithin(effs_combined; f=(x, y) -> -(x, y))
    within_combined = @chain effs_combined begin
        unstack(:sampler, :eff)
        @rtransform(
            :ΔUS_SR = f(:UncertaintySampling, :SimpleRandom),
            :ΔUS_WBA = f(:UncertaintySampling, :WeightedBalancedAcceptance),
            :ΔWBA_SR = f(:WeightedBalancedAcceptance, :SimpleRandom),
            :ΔRI_SR = f($("Realized interactions"), $("Species richness")),
            :ΔRI_FR = f($("Realized interactions"), $("Focal species range")),
            :ΔFR_SR = f($("Focal species range"), $("Species richness")),
        )
        select(:sim, :set, :occ, r"Δ")
        stack(r"Δ")
        dropmissing()
    end
    return within_combined
end
within_combined = comparewithin(effs_combined)
within_combined_ndi = comparewithin(effs_combined; f=ndi)
within_combined_log = comparewithin(effs_combined; f=(x, y) -> log(x) - log(y))

# Visualize
begin
    m = mapping(:variable, :value => "Δefficiency"; color=:value => (x -> x >= 0.0))
    rains = visual(
        RainClouds;
        markersize=7,
        jitter_width=0.15,
        plot_boxplots=false,
        clouds=nothing,
        orientation=:horizontal,
    )
    vline = mapping([0]) * visual(VLines; linestyle=:dash)
end
let d = within_combined_log
    Random.seed!(42)
    d1 = @rsubset(d, :set == "Samplers")
    d2 = @rsubset(d, :set == "Layers")
    m = mapping(:variable, :value => "efficiency (log)"; color=:value => (x -> x >= 0.0))
    f = Figure()
    fg1 = draw!(f[1, 1], data(d1) * m * rains + vline; axis=(; title="Samplers"))
    fg2 = draw!(f[2, 1], data(d2) * m * rains + vline; axis=(; title="Layers"))
    linkxaxes!(fg1..., fg2...)
    f
end
save(plotsdir("efficiency_comparison_log.png"), current_figure())

@chain begin
    groupby(within_combined, [:set, :variable])
    @combine(:prop = sum((:value .>= 0) ./ length(:value)))
end

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

# Rename
effs_samplers.sampler =
    replace.(
        effs_samplers.sampler,
        "UncertaintySampling" => "Uncertainty Sampling",
        "WeightedBalancedAcceptance" => "Weighted Balanced Acceptance",
        "SimpleRandom" => "Simple Random",
    )

# Combine for convenience
effs_combined = vcat(effs_samplers, effs_optimized)

# Define color sets
cols = [
    # Interaction types
    "possible" => Makie.wong_colors()[2],
    "realized" => Makie.wong_colors()[3],
    "detected" => Makie.wong_colors()[4],
    # Samplers
    "Uncertainty Sampling" => Makie.wong_colors()[2],
    "Weighted Balanced Acceptance" => Makie.wong_colors()[3],
    "Simple Random" => Makie.wong_colors()[1],
    # Layers
    "Focal species range" => Makie.wong_colors()[2],
    "Species richness" => Makie.wong_colors()[4],
    "Realized interactions" => Makie.wong_colors()[5],
    "Probabilistic range" => Makie.wong_colors()[6],
]

## Efficiency only

# Define sorting order across figures
sortedsamplers = ["Uncertainty Sampling", "Weighted Balanced Acceptance", "Simple Random"]
sortedlayers = [
    "Focal species range",
    "Realized interactions",
    "Probabilistic range",
    "Species richness",
]
sortedlayout = [sortedsamplers..., sortedlayers...]

# Violin
begin
    Random.seed!(42) # for jitter
    eff = :eff => "efficiency"
    layer =
        mapping(:sampler => sorter(sortedlayout) => "", eff; color=:sampler) *
        visual(RainClouds; markersize=5, jitter_width=0.1, plot_boxplots=false)
    f1 = data(effs_samplers) * layer
    f2 = data(effs_optimized) * layer
    f = Figure(; size=(700, 700))
    sc = scales(; Color=(; palette=cols))
    ylog2f = (; ytickformat=vs -> [rich("2", superscript("$(Int(v))")) for v in vs])
    ylog2 = (; axis=(; yticks=2:2:100))
    fg1 = draw!(f[1, 1], f1, sc; ylog2...)
    fg2 = draw!(f[2, 1], f2, sc; ylog2...)
    linkyaxes!(fg1..., fg2...)
    pad = (-45, 0, 20, 0)
    Label(f[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(f[2, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)
    f
end
save(plotsdir("efficiency_distribution.png"), f)

# Efficiency given species degree
begin
    Random.seed!(42)
    ranknames = renamer(1 => "1.0", 2 => "0.66", 3 => "0.33", 4 => "0.07")
    ranks = :rank => ranknames => "Within-simulation Degree Percentile Rank"
    sppalette = from_continuous(cgrad(:viridis; rev=true))
    f = Figure(; size=(650, 650))
    fg1 = draw!(
        f[1, 1:2],
        data(effs_species) *
        mapping(ranks, eff; color=ranks) *
        visual(RainClouds; markersize=5, jitter_width=0.1, plot_boxplots=false),
        scales(; Color=(; palette=sppalette));
        axis=(; yticks=2:4:100),
    )
    fg2 = draw!(
        f[2, 1],
        data(effs_species) *
        mapping(:deg => "Degree", eff; color=:rank => ranknames => "Percentile rank") *
        visual(Scatter),
        scales(; Color=(; palette=sppalette));
        ylog2...,
    )
    legend!(f[2, 2], fg2)
    f
end
save(plotsdir("efficiency_distribution_species.png"), f)

## Efficiency & Occupancy

# Scatter & smooth
begin
    occ = :occ => "occupancy"
    layout = mapping(occ, eff; color=:sampler) * (visual(Scatter) + linear())
    scale = scales(; Color=(; palette=cols))
    legend = (; position=:bottom)
    f1 = data(effs_samplers) * layout * mapping(; col=:sampler => sorter(sortedlayout))
    f2 = data(effs_optimized) * layout * mapping(; col=:sampler => sorter(sortedlayout))
end
draw(f1, scale; legend=legend)
draw(f2, scale; legend=legend)
fig = draw(
    data(effs_combined) * layout * mapping(; layout=:sampler => sorter(sortedlayout)),
    scales(; Color=(; palette=cols, legend=false));
    figure=(; size=(800, 450)),
)

# Same with independent subfigures
fig = let
    f = Figure(; size=(800, 500))
    fg1 = draw!(f[1, 1:3], f1, scale; ylog2...)
    fg2 = draw!(f[2, 1:4], f2, scale; ylog2...)
    linkaxes!(fg1..., fg2...)
    pad = (-50, 0, 30, 0)
    Label(f[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(f[2, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)
    f
end
save(plotsdir("efficiency_occupancy.png"), fig)

# Group columns in single panel
fig = let
    f = Figure(; size=(800, 450))
    fg1 = draw!(f[1, 1], f1 * mapping(; col=:set), scale; ylog2...)
    legend!(f[2, 1], fg1; position=:bottom, tellheight=true, tellwidth=false)
    fg2 = draw!(f[1, 2], f2 * mapping(; col=:set), scale; ylog2...)
    legend!(f[2, 2], fg2; position=:bottom, tellheight=true, tellwidth=false)
    linkyaxes!(fg1..., fg2...)
    f
end
save(plotsdir("_xtras/", "efficiency_occupancy_scatter.png"), fig)

## Efficiency & Occupancy with species degree

# Scatter & smooth with percentile rank
sppalette = from_continuous(cgrad(:viridis; rev=true))
fig = draw(
    data(effs_species) * layout * mapping(; color=ranks, col=ranks),
    scales(; Color=(; palette=sppalette));
    figure=(; size=(800, 300)),
    legend=legend,
    axis=(; yticks=4:4:100),
)
save(plotsdir("efficiency_occupancy_species.png"), fig)

# Group all ranks in single figure
fig = draw(
    data(effs_species) * layout * mapping(; color=ranks),
    scales(; Color=(; palette=sppalette));
    legend=legend,
    ylog2...,
)
save(plotsdir("_xtras", "efficiency_occupancy_species_rank.png"), fig)

# Species degree & efficiency-occupancy
fig = draw(data(effs_species) * layout * mapping(; color=:deg); legend=legend, ylog2...)
save(plotsdir("_xtras/", "efficiency_occupancy_species_degree.png"), fig)

## Within-simulation comparison

# NDI
ndi(x, y) = (x - y) / (x + y)

# Efficiency difference
function efficiency_difference(n, n2; k=10_000)
    n = exp(n)
    n2 = exp(n2)
    return ((n * log(n) - n * log(n + k) + k) - (n2 * log(n2) - n2 * log(n2 + k) + k))
end

# Separate results per simulation
function comparewithin(effs_combined; f=(x, y) -> -(x, y))
    within_combined = @chain effs_combined begin
        unstack(:sampler, :eff)
        @rtransform(
            :ΔUS_SR = f($("Uncertainty Sampling"), $("Simple Random")),
            :ΔUS_WBA = f($("Uncertainty Sampling"), $("Weighted Balanced Acceptance")),
            :ΔWBA_SR = f($("Weighted Balanced Acceptance"), $("Simple Random")),
            :ΔRI_SR = f($("Realized interactions"), $("Species richness")),
            :ΔRI_FR = f($("Realized interactions"), $("Focal species range")),
            :ΔFR_SR = f($("Focal species range"), $("Species richness")),
            :ΔRI_PR = f($("Realized interactions"), $("Probabilistic range")),
            :ΔFR_PR = f($("Focal species range"), $("Probabilistic range")),
            :ΔPR_SR = f($("Probabilistic range"), $("Species richness")),
        )
        select(:sim, :set, :occ, r"Δ")
        stack(r"Δ")
        dropmissing()
    end
    return within_combined
end
within_combined = comparewithin(effs_combined)
within_combined_ndi = comparewithin(effs_combined; f=ndi)
within_combined_log = comparewithin(effs_combined; f=(x, y) -> (x / y))
within_combined_dif = comparewithin(
    effs_combined; f=(n, n2) -> efficiency_difference(n, n2; k=10_000)
)

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
    vline = mapping([0.5]) * visual(VLines; linestyle=:dash)
end
let d = within_combined_dif
    Random.seed!(42)
    d1 = @rsubset(d, :set == "Samplers")
    d2 = @rsubset(d, :set == "Layers")
    m = mapping(
        :variable => "comparison",
        :value => "Efficiency difference";
        color=:value => (x -> x <= 0.0),
    )
    f = Figure()
    xlog2f = vs -> [rich("2", superscript("$(v)")) for v in vs]
    xlog2 = (; axis=(; xtickformat=xlog2f))
    # xaxis = (; axis=(; xticks=0:2:10))
    fg1 = draw!(f[1, 1], data(d1) * m * rains + vline;)
    fg2 = draw!(f[2:3, 1], data(d2) * m * rains + vline;)
    linkxaxes!(fg1..., fg2...)
    pad = (-100, 0, 10, 0)
    Label(f[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(f[2, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)
    f
end
save(plotsdir("efficiency_comparison.png"), current_figure())

# Reduce number of comparisons
logit(p) = 1 / (1 + exp(-p))
comps_dict = Dict(
    "ΔUS_SR" => "Simple Random",
    "ΔUS_WBA" => "Weighted Balanced Acceptance",
    "ΔRI_FR" => "Realized Interactions",
    "ΔFR_SR" => "Species Richness",
    "ΔFR_PR" => "Probabilistic Range",
)
within_combined_dif2 = @chain within_combined_dif begin
    @rsubset(:variable in keys(comps_dict))
    @rtransform(:variable = comps_dict[:variable])
    @rtransform(:value = :variable == "Realized Interactions" ? -(:value) : :value)
    @rtransform(:value = -:value)
    # @rtransform(:value = logit(:value))
end
let d = within_combined_dif2
    Random.seed!(42)
    d1 = @rsubset(d, :set == "Samplers")
    d2 = @rsubset(d, :set == "Layers")
    m = mapping(
        :variable => "",
        :value => logit => "Efficiency compared to reference (Uncertainty Sampling)";
        color=:value => (x -> x >= 0.0),
    )
    f = Figure()
    xlog2f = vs -> [rich("2", superscript("$(v)")) for v in vs]
    xlog2 = (; axis=(; xtickformat=xlog2f))
    # xaxis = (; axis=(; xticks=0:2:10))
    fg1 = draw!(f[1, 1], data(d1) * m * rains + vline;)
    fg2 = draw!(
        f[2:3, 1],
        data(d2) * m * rains + vline;
        axis=(; xlabel="Efficiency compared to reference (Focal Range)"),
    )
    linkxaxes!(fg1..., fg2...)
    pad = (-100, 0, 10, 0)
    Label(f[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(f[2, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)
    f
end
save(plotsdir("efficiency_comparison_reduced.png"), current_figure())

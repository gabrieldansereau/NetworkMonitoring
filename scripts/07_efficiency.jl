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

## Summarize efficiency simulations

# Summmarize results not combined previously
function summarize_monitored(df)
    monitored = @chain df begin
        groupby([:sp, :type, :sampler, :nbon])
        @combine(
            :low = quantile(:monitored, 0.05),
            :med = median(:monitored),
            :upp = quantile(:monitored, 0.95),
            :deg = maximum(:deg)
        )
        @rtransform(:low = :low / :deg, :med = :med / :deg, :upp = :upp / :deg,)
    end
    return monitored
end

# Use job id to vary parameters
files = filter(startswith("monitored_samplers"), readdir(datadir("efficiency")))
ids = sort(parse.(Int, replace.(files, "monitored_samplers-" => "", ".csv" => "")))
sims_samplers = DataFrame()
sims_optimized = DataFrame()
sims_species = DataFrame()
for id in ids
    # Load all results
    idp = lpad(id, 2, "0")
    monitored_samplers_all = CSV.read(
        datadir("efficiency", "monitored_samplers-$idp.csv"), DataFrame
    )
    monitored_optimized_all = CSV.read(
        datadir("efficiency", "monitored_optimized-$idp.csv"), DataFrame
    )
    monitored_species_all = CSV.read(
        datadir("efficiency", "monitored_spp-$idp.csv"), DataFrame
    )

    # Summmarize results not combined previously
    monitored_samplers = summarize_monitored(monitored_samplers_all)
    monitored_optimized = summarize_monitored(monitored_optimized_all)
    monitored_species = summarize_monitored(monitored_species_all)

    # Add species rank
    monitored_species_occ = CSV.read(
        datadir("efficiency", "monitored_spp_occ-$idp.csv"), DataFrame
    )
    leftjoin!(monitored_species, monitored_species_occ; on=:sp)

    # Add sim id
    @select!(monitored_samplers, :sim = id, All())
    @select!(monitored_optimized, :sim = id, All())
    @select!(monitored_species, :sim = id, All())

    # Collect
    append!(sims_samplers, monitored_samplers)
    append!(sims_optimized, monitored_optimized)
    append!(sims_species, monitored_species)
end
sims_samplers
sims_optimized
sims_species

# Export
CSV.write(datadir("sims_efficiency_samplers.csv"), sims_samplers)
CSV.write(datadir("sims_efficiency_optimized.csv"), sims_optimized)
CSV.write(datadir("sims_efficiency_species.csv"), sims_species)

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

# Select UncertaintySampling only as example
sims_set = @rsubset(sims_samplers, :sampler == "UncertaintySampling")
@rsubset!(sims_set, :sim <= 10)

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

## Occupancy

# Get layer occupancy
occupancy(l) = length(findall(isone, l)) / length(l)

# Calculate occupancy
ids = sort(unique(sims_samplers.sim))
occup = zeros(length(ids))
for i in ids
    idp = lpad(i, 2, "0")
    l = SDT.SDMLayer(datadir("efficiency", "layer_sp_range-$idp.tiff"))
    occup[i] = occupancy(l)
end
occupdf = DataFrame(; sim=ids, occ=occup)

# Calculate efficiency & assign occupancy
effs_samplers = @chain sims_samplers begin
    @groupby(:sim, :sampler)
    @combine(:eff = efficiency(:nbon, :med))
    leftjoin(occupdf; on=:sim)
    @transform(:occ = occup[:sim], :set = "Samplers")
end
effs_optimized = @chain sims_optimized begin
    @groupby(:sim, :sampler)
    @combine(:eff = efficiency(:nbon, :med))
    leftjoin(occupdf; on=:sim)
    @transform(:set = "Layers")
end
effs_species = @chain sims_species begin
    @groupby(:sim, :sampler, :sp, :deg, :rank, :occ)
    @combine(:eff = efficiency(:nbon, :med))
    @transform(:set = "Species")
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
sortedsamplers = ["UncertaintySampling", "WeightedBalancedAcceptance", "SimpleRandom"]
sortedlayers = ["Focal species range", "Species richness", "Realized interactions"]
sortedlayout = [sortedsamplers..., sortedlayers...]
begin
    occ = :occ => "occupancy"
    efflog = :eff => log => "log(efficiency)"
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
save(plotsdir("saturation_occupancy_scatter_facets.png"), fig)

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
save(plotsdir("saturation_occupancy_scatter.png"), fig)

# Species degree & occupancy
fig = draw(data(effs_species) * layout * mapping(; color=:deg); legend=legend)
save(plotsdir("saturation_occupancy_degree.png"), fig)

# Species rank instead of degree
ranknames = renamer(1 => "1.0", 2 => "0.66", 3 => "0.33", 4 => "0.07")
fig = draw(
    data(effs_species) *
    layout *
    mapping(; color=:rank => ranknames => "Within-simulation Degree Percentile Rank"),
    scales(; Color=(; palette=from_continuous(cgrad(:viridis; rev=true))));
    legend=legend,
)
save(plotsdir("saturation_occupancy_degree_ranked.png"), fig)

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
        visual(RainClouds; markersize=5, jitter_width=0.1, plot_boxplots=false)
    f1 = data(effs_samplers) * layer
    f2 = data(effs_optimized) * layer
    f = Figure(; size=(700, 700))
    draw!(f[1, 1], f1, scales(; Color=(; palette=cols)))
    draw!(f[2, 1], f2, scales(; Color=(; palette=cols)))
    f
end
save(plotsdir("saturation_efficiency_distribution.png"), f)

## Within-simulation comparison

# Separate results per simulation
within_samplers = @chain effs_samplers begin
    unstack(:sampler, :eff)
    @rtransform(
        :ΔUS_SR = :UncertaintySampling - :SimpleRandom,
        :ΔUS_WBA = :UncertaintySampling - :WeightedBalancedAcceptance,
        :ΔWBA_SR = :WeightedBalancedAcceptance - :SimpleRandom,
    )
    select(:sim, :set, :occ, r"Δ")
    stack(Not([:sim, :set, :occ]))
end
within_samplers = @chain effs_optimized begin
    unstack(:sampler, :eff)
    @rtransform(
        :ΔSR_RI = $("Species richness") - $("Realized interactions"),
        :ΔUS_RI = $("Focal species range") - $("Realized interactions"),
        :ΔUS_SR = $("Focal species range") - $("Species richness"),
    )
    select(:sim, :set, :occ, r"Δ")
    stack(Not([:sim, :set, :occ]))
end

# Visualize
begin
    fig =
        data(within_samplers) *
        # mapping(:variable, :value => "Δefficiency"; color=:occ) *
        mapping(:variable, :value => "Δefficiency"; color=:value => (x -> x >= 0.0)) *
        visual(
            RainClouds; markersize=7, jitter_width=0.15, plot_boxplots=false, clouds=nothing
        )
    draw(fig + data((; y=[0])) * mapping(:y) * visual(HLines; linestyle=:dash))
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

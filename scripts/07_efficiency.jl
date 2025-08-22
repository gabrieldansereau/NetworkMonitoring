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
function efficiency(x, y; A=LinRange(-5.0, 15.0, 10_000))
    err = zeros(length(A))
    for i in eachindex(A)
        f = saturation(exp2(A[i]))
        err[i] = sqrt(sum((y .- f(x)) .^ 2.0))
    end
    return exp2(A[last(findmin(err))])
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

# Export
CSV.write(datadir("efficiency_samplers.csv"), effs_samplers)
CSV.write(datadir("efficiency_optimized.csv"), effs_optimized)
CSV.write(datadir("efficiency_species.csv"), effs_species)

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
save(plotsdir("efficiency_within_example.png"), fig)
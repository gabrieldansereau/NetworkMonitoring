# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module

## Summarize efficiency simulations

# Use job id to vary parameters
files = filter(startswith("monitored_samplers"), readdir(datadir("efficiency")))
ids = sort(parse.(Int, replace.(files, "monitored_samplers-" => "", ".csv" => "")))
sets = ["samplers", "optimized", "spp", "estimations"]
sims_dict = Dict()
for set in sets
    # Create data frame for set results
    sims_dict[set] = DataFrame()
    for id in ids
        # Load all results
        idp = lpad(id, 2, "0")
        file = datadir("efficiency", "monitored_$set-$idp.csv")
        isfile(file) || continue
        monitored_all = CSV.read(datadir(file), DataFrame)

        # Summmarize results not combined previously
        monitored = summarize_focal(monitored_all; id=id)

        # Add species rank and occupancy
        if set == "spp"
            monitored_occ = CSV.read(
                datadir("efficiency", "monitored_spp_occ-$idp.csv"), DataFrame
            )
            leftjoin!(monitored, monitored_occ; on=[:sp, :sim])
        end

        # Collect
        append!(sims_dict[set], monitored)
    end
end
sims_samplers = sims_dict["samplers"]
sims_optimized = sims_dict["optimized"]
sims_species = sims_dict["spp"]
sims_estimations = sims_dict["estimations"]

## Efficiency

# Select UncertaintySampling only as example
sims_set = @rsubset(sims_samplers, :sampler == "Uncertainty Sampling")
@rsubset!(sims_set, :sim <= 10)

# Visualize
begin
    f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel="Number of sites",
        ylabel="Proportion of sampled interactions",
        title="Illustration of replicates & efficiency - Uncertainty Sampling",
    )
    for iter in unique(sims_set.sim)
        X₁ = @rsubset(sims_set, :sim == iter)
        x = X₁.nbon
        y = X₁.med
        eff = efficiency_gridsearch(x, y)
        scatter!(ax, x, y; color=:lightgrey)
        lines!(ax, x, saturation(eff)(x); color=:black, linestyle=:dash)
    end
    f
end
save(plotsdir("efficiency_example.png"), f)

## Occupancy

# Calculate occupancy
ids = sort(unique(sims_samplers.sim))
occup = Dict(ids .=> zeros(length(ids)))
for i in ids
    idp = lpad(i, 2, "0")
    l = SDT.SDMLayer(datadir("efficiency", "layer_sp_range-$idp.tiff"))
    occup[i] = occupancy(l)
end
occupdf = DataFrame(; sim=ids, occ=[occup[i] for i in ids])

# Calculate efficiency & assign occupancy
effs_species = @chain sims_species begin
    @groupby(:sim, :set, :sp, :deg, :rank, :occ)
    @combine(:eff = efficiency(:nbon, :med))
    rename(:sp => :variable)
end
effs_samplers = @chain sims_samplers begin
    @groupby(:sim, :set, :sampler)
    @combine(:eff = efficiency(:nbon, :med))
    @rtransform(:occ = occup[:sim])
    rename(:sampler => :variable)
end
effs_optimized = @chain sims_optimized begin
    @groupby(:sim, :set, :layer)
    @combine(:eff = efficiency(:nbon, :med))
    @rtransform(:occ = occup[:sim])
    rename(:layer => :variable)
end
effs_estimations = @chain sims_estimations begin
    @groupby(:sim, :set, :layer)
    @combine(:eff = efficiency(:nbon, :med))
    @rtransform(:occ = occup[:sim])
    rename(:layer => :variable)
end

# Export
CSV.write(datadir("efficiency_samplers.csv"), effs_samplers);
CSV.write(datadir("efficiency_optimized.csv"), effs_optimized);
CSV.write(datadir("efficiency_species.csv"), effs_species);
CSV.write(datadir("efficiency_estimations.csv"), effs_estimations);

## Within-simulation variation

# Load one simulation example
monitored_samplers = CSV.read(datadir("monitored_samplers.csv"), DataFrame)
@rsubset!(monitored_samplers, :sampler == "Uncertainty Sampling")

# Visualize
let b = monitored_samplers
    eff = efficiency_gridsearch(b.nbon, b.med)
    band(b.nbon, b.low, b.upp; alpha=0.4, color=Makie.wong_colors()[2], label="intervals")
    lines!(b.nbon, b.med; color=Makie.wong_colors()[2], linewidth=1.5, label="median")
    lines!(
        b.nbon, saturation(eff)(b.nbon); color=:black, linestyle=:dash, label="efficiency"
    )
    hlines!([1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb)
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
        eff = efficiency_gridsearch(bi.nbon, bi.monitored; A=LinRange(-12.0, 12.0, 5000))
        lines!(
            bi.nbon,
            saturation(eff)(bi.nbon);
            color=:grey,
            linestyle=:solid,
            label="individual saturation curves",
        )
    end
    eff = efficiency_gridsearch(b.nbon, b.med)
    band!(
        b.nbon,
        b.low,
        b.upp;
        alpha=0.4,
        color=colours["Uncertainty Sampling"],
        label="90% percentile",
    )
    scatter!(b.nbon, b.med; color=colours["Uncertainty Sampling"], label="median points")
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
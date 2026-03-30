# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module

## Summarize efficiency simulations

# Use job id to vary parameters
sets = ["samplers", "optimized", "spp", "estimations"]
sims_dict = Dict()
for set in sets
    # List available results
    files = filter(startswith("monitored_$set"), readdir(datadir("efficiency")))
    ids = sort(
        replace.(files, "monitored_$set-" => "", ".csv" => ""); by=x -> parse(Int, x)
    )
    # Create data frame for set results
    sims_dict[set] = DataFrame()
    for id in ids
        # Load all results
        file = datadir("efficiency", "monitored_$set-$id.csv")
        isfile(file) || continue
        monitored_all = CSV.read(datadir(file), DataFrame)

        # Separate results with missing values
        monitored_missings = filter(:monitored => ismissing, monitored_all)
        filter!(:monitored => !ismissing, monitored_all)

        # Summmarize results not combined previously
        _confint = set == "estimations" ? true : false
        monitored = summarize_focal(
            monitored_all; id=parse(Int, id), confint=_confint, α=0.10
        )

        # Add species rank and occupancy
        if set == "spp"
            monitored_occ = CSV.read(
                datadir("efficiency", "monitored_sp_occ-$id.csv"), DataFrame
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
        eff = efficiency_gridsearch(x, y; f=exp)
        scatter!(ax, x, y; color=:lightgrey)
        lines!(ax, x, saturation(eff)(x); color=:black, linestyle=:dash)
    end
    f
end
save(plotsdir("supp", "efficiency_example.png"), f)

## Occupancy

# Calculate occupancy
files = filter(startswith("layer_sp_range"), readdir(datadir("efficiency")))
sort!(files; by=x -> parse(Int, replace(x, "layer_sp_range-" => "", ".tiff" => "")))
occup = Dict()
for (i, f) in enumerate(files)
    l = SDT.SDMLayer(datadir("efficiency", f); bandnumber=1)
    occup[i] = occupancy(l)
end
occupdf = DataFrame(; sim=collect(keys(occup)), occ=collect(values(occup)))
sort!(occupdf, :sim)

# Calculate efficiency & assign occupancy
effs_species = @chain sims_species begin
    @groupby(:sim, :set, :sp, :deg, :rank, :occ)
    @combine(:eff = efficiency(:nbon, :med; f=exp))
    rename(:sp => :variable)
end
effs_samplers = @chain sims_samplers begin
    @groupby(:sim, :set, :sampler)
    @combine(:eff = efficiency(:nbon, :med; f=exp))
    @rtransform(:occ = occup[:sim])
    rename(:sampler => :variable)
end
effs_optimized = @chain sims_optimized begin
    @groupby(:sim, :set, :layer)
    @combine(:eff = efficiency(:nbon, :med; f=exp))
    @rtransform(:occ = occup[:sim])
    rename(:layer => :variable)
end
effs_estimations = @chain sims_estimations begin
    @groupby(:sim, :set, :layer)
    @combine(
        :eff = efficiency(:nbon, :med; f=exp),
        :eff_low = efficiency(:nbon, :confint_low; f=exp),
        :eff_upp = efficiency(:nbon, :confint_upp; f=exp)
    )
    @rtransform(:occ = occup[:sim])
    rename(:layer => :variable)
end

# Rename estimations with same number of digits
@transform!(
    effs_estimations,
    :variable = [
        @sprintf("%s-%.2f", s[1], parse(Float64, s[2])) for s in split.(:variable, "-")
    ]
)

# Add the offset as a column
@chain effs_estimations begin
    @rtransform!(
        :offset = parse(
            Float64, replace(:variable, "Over-" => "", "Under" => "", "True-" => "")
        )
    )
    select!(:sim, :set, :variable, :offset, All())
    sort!(:sim)
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

# Load pre-summary results
file1 = filter(startswith("monitored_samplers"), readdir(datadir("efficiency")))[1]
sim1_samplers = @chain begin
    CSV.read(datadir("efficiency", file1), DataFrame)
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
        eff = efficiency_gridsearch(
            bi.nbon, bi.monitored; A=LinRange(-12.0, 12.0, 5000), f=exp
        )
        lines!(
            bi.nbon,
            saturation(eff)(bi.nbon);
            color=:grey,
            linestyle=:solid,
            label="individual saturation curves",
        )
    end
    eff = efficiency_gridsearch(b.nbon, b.med; f=exp)
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
save(plotsdir("supp", "efficiency_within_example.png"), fig)
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
        monitored = summarize_focal(monitored_all; id=parse(Int, id), confint=true, α=0.10)

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

## Efficiency example

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

## Efficiency for species, samplers and optimized simulations

# Calculate efficiency & assign occupancy
eff_opt = (; f=exp, option=:n_at_p, p=0.8)
gp_vars = [:sim, :set]
ordered_vars = [:sim, :variable, :eff, :eff_low, :eff_upp, :rmse, :occ, :set]
effs_species = @chain sims_species begin
    @groupby(gp_vars, :sp, :deg, :rank, :occ)
    @combine(
        :eff = Ref(efficiency(:nbon, :med; eff_opt..., rmse=true)),
        :eff_low = efficiency(:nbon, :confint_low; eff_opt...),
        :eff_upp = efficiency(:nbon, :confint_upp; eff_opt...),
    )
    @rtransform(:rmse = getfield(:eff, ^(:rmse)), :eff = getfield(:eff, ^(:ei)))
    rename(:sp => :variable)
    select(ordered_vars, All())
end
effs_samplers = @chain sims_samplers begin
    @groupby(gp_vars, :sampler)
    @combine(
        :eff = Ref(efficiency(:nbon, :med; eff_opt..., rmse=true)),
        :eff_low = efficiency(:nbon, :confint_low; eff_opt...),
        :eff_upp = efficiency(:nbon, :confint_upp; eff_opt...),
    )
    @rtransform(:rmse = getfield(:eff, ^(:rmse)), :eff = getfield(:eff, ^(:ei)))
    @rtransform(:occ = occup[:sim])
    rename(:sampler => :variable)
    select(ordered_vars, All())
end
effs_optimized = @chain sims_optimized begin
    @groupby(gp_vars, :layer)
    @combine(
        :eff = Ref(efficiency(:nbon, :med; eff_opt..., rmse=true)),
        :eff_low = efficiency(:nbon, :confint_low; eff_opt...),
        :eff_upp = efficiency(:nbon, :confint_upp; eff_opt...),
    )
    @rtransform(:rmse = getfield(:eff, ^(:rmse)), :eff = getfield(:eff, ^(:ei)))
    @rtransform(:occ = occup[:sim])
    rename(:layer => :variable)
    select(ordered_vars, All())
end

# Export
CSV.write(datadir("efficiency_samplers.csv"), effs_samplers);
CSV.write(datadir("efficiency_optimized.csv"), effs_optimized);
CSV.write(datadir("efficiency_species.csv"), effs_species);

## Efficiency for range estimations

# Prior transformations and grouping for range estimations
gdf = @chain sims_estimations begin
    @rtransform(:degmax = :degmax / :deg)
    rename(:degmax => :pmax)
    @groupby(:sim, :layer, :offset)
end

# Calculate the efficiency for the range estimations
effs_estimations = DataFrame()
pmax_opt = :n_at_pmax4
@showprogress "Efficiency for range estimations" for gd in gdf
    # Extract group info
    sim = first(gd.sim)
    layer = first(gd.layer)
    offset = first(gd.offset)
    deg = first(gd.deg)
    pmax = first(gd.pmax)
    occ = occup[sim]
    set = first(gd.set)
    # Compute efficiency
    eff, rmse = efficiency(
        gd.nbon, gd.med; f=exp, pmax=pmax, p=0.80, option=pmax_opt, rmse=true
    )
    eff_low = efficiency(gd.nbon, gd.confint_low; f=exp, pmax=pmax, p=0.80, option=pmax_opt)
    eff_upp = efficiency(gd.nbon, gd.confint_upp; f=exp, pmax=pmax, p=0.80, option=pmax_opt)
    # Export
    row = (;
        sim=sim,
        variable=layer,
        offset=offset,
        eff=eff,
        eff_low=eff_low,
        eff_upp=eff_upp,
        rmse=rmse,
        deg=deg,
        pmax=pmax,
        occ=occ,
        set=set,
    )
    push!(effs_estimations, row)
end
effs_estimations

# Export
CSV.write(datadir("efficiency_estimations.csv"), effs_estimations);

## Adjust efficiency estimation for point density

# Reference values
nbon_max = maximum(sims_estimations.nbon)
nbon_ref = maximum(@rsubset(sims_estimations, :offset == 0.0).nbon)
ref_file = first(filter(startswith("layer_sp_range"), readdir(datadir("efficiency"))))
ntot = length(SDT.SDMLayer(datadir("efficiency", ref_file); bandnumber=1))
ref_density = Dict(k => nbon_ref / (ntot * v) for (k, v) in occup)

# Adjust sampling effort
sims_effort = @chain sims_estimations begin
    # Remove unnecessary extreme simulations
    @rsubset(:offset <= 0.5)
    # Prepare the columns
    @select(:sim, :layer, :offset, :nbon)
    @rtransform(:occ = occup[:sim])
    # Measure the total area given the offset
    @rtransform(:area = (1 + :offset) * :occ)
    @rtransform(:area = :area > 1.0 ? 1.0 : :area)
    # Get the sampling density and effort-adjusted number of sites
    @rtransform(
        :sampling_density = :nbon / (ntot * :area),
        :ref_density = ref_density[:sim],
        :nmax = round(Int, nbon_ref * (1 + :offset))
    )
    @rtransform(:neff = :nbon * nbon_max / :nmax)
    # Arrange columns
    @rselect(:sim, :layer, :offset, :nbon, :nmax, :neff, Not(:area))
    # Select the sampling-adjusted observations
    # @rsubset(:sampling_density <= :ref_density)
    @rsubset(:neff <= nbon_max)
    # Join back with monitoring resuts
    leftjoin(sims_estimations; on=[:sim, :layer, :offset, :nbon])
end

# Calculate efficiency on effort-adjusted n
effs_effort = @chain sims_effort begin
    @groupby(:sim, :set, :layer, :offset)
    @combine(
        :eff = efficiency(:nbon, :med; f=exp),
        :eff_low = efficiency(:nbon, :confint_low; f=exp),
        :eff_upp = efficiency(:nbon, :confint_upp; f=exp)
    )
    @rtransform(:occ = occup[:sim])
    rename(:layer => :variable)
end

# Export
CSV.write(datadir("efficiency_estimations-effort_adjusted.csv"), effs_effort);

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
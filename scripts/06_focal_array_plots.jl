using DrWatson
@quickactivate :NetworkMonitoring

update_theme!(; CairoMakie=(; px_per_unit=2.0))

## Load focal species results

# Use job id to vary parameters
id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
idp = lpad(id, 2, "0")

# Load all results
monitored_samplers_all = CSV.read(datadir("monitored_samplers-$idp.csv"), DataFrame)
monitored_optimized_all = CSV.read(datadir("monitored_optimized-$idp.csv"), DataFrame)

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
        rename(:type => :var)
    end
    return monitored
end
monitored_samplers = summarize_monitored(monitored_samplers_all)
monitored_optimized = summarize_monitored(monitored_optimized_all)

# Get species with highest degree
sp = monitored_samplers.sp[1]
deg = maximum(monitored_samplers.deg)

# Load layers used for optimization
focal_sp_range = SDT.SDMLayer(datadir("layer_sp_range-$idp.tiff"))
richness_spp = SDT.SDMLayer(datadir("layer_richness_spp-$idp.tiff"))
degree_realized = SDT.SDMLayer(datadir("layer_degree_realized-$idp.tiff"))

## Load all simulations

# Use job id to vary parameters
ids = 1:10
sims_samplers = DataFrame()
sims_optimized = DataFrame()
for id in ids
    # Load all results
    idp = lpad(id, 2, "0")
    monitored_samplers_all = CSV.read(datadir("monitored_samplers-$idp.csv"), DataFrame)
    monitored_optimized_all = CSV.read(datadir("monitored_optimized-$idp.csv"), DataFrame)

    # Summmarize results not combined previously
    monitored_samplers = summarize_monitored(monitored_samplers_all)
    monitored_optimized = summarize_monitored(monitored_optimized_all)

    # Add sim id
    @select!(monitored_samplers, :sim = idp, All())
    @select!(monitored_optimized, :sim = idp, All())

    # Collect
    append!(sims_samplers, monitored_samplers)
    append!(sims_optimized, monitored_optimized)
end
sims_samplers
sims_optimized

# Re-summarize
monitored_samplers = @chain sims_samplers begin
    @rtransform(:lowdiff = :med - :low, :uppdiff = :upp - :med)
    groupby([:sampler, :nbon])
    @combine(
        :low = mean(:lowdiff),
        :med = median(:med),
        :upp = mean(:uppdiff),
        :lowmed = median(:lowdiff),
        :uppmed = median(:uppdiff),
    )
end
monitored_optimized = @chain sims_optimized begin
    @rtransform(:lowdiff = :med - :low, :uppdiff = :upp - :med)
    groupby([:sampler, :nbon])
    @combine(
        :low = mean(:lowdiff),
        :med = median(:med),
        :upp = mean(:uppdiff),
        :lowmed = median(:lowdiff),
        :uppmed = median(:uppdiff),
    )
end

## Define labels and colors for all plots

# Colors
cols = Dict{Any,Any}(
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
)

## Explore variations with different sampler

# Generate BON examples
begin
    Random.seed!(42)
    bons = Dict()
    bons["UncertaintySampling"] = BON.sample(BON.UncertaintySampling(100), focal_sp_range)
    bons["WeightedBalancedAcceptance"] = BON.sample(
        BON.WeightedBalancedAcceptance(100), focal_sp_range
    )
    bons["SimpleRandom"] = BON.sample(BON.SimpleRandom(100), focal_sp_range)
end

# Plot
begin
    res = filter(:sampler => ==("UncertaintySampling"), sims_samplers)
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:100:500
    )
    ax1 = Axis(
        gb[1, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax2 = Axis(
        gb[2, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax3 = Axis(
        gb[3, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        for sim in unique(res.sim)
            bsim = filter(:sim => ==(sim), b)
            band!(ax, bsim.nbon, bsim.low, bsim.upp; alpha=0.3, label=s, color=cols[s])
            lines!(ax, bsim.nbon, bsim.med; label=s, color=cols[s])
        end
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Legend(ga[2,1], ax, orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, focal_sp_range)
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_array_samplers_one.png"), fig)

# Plot
begin
    res = sims_samplers
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:100:500
    )
    ax1 = Axis(
        gb[1, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax2 = Axis(
        gb[2, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax3 = Axis(
        gb[3, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        for sim in unique(res.sim)
            bsim = filter(:sim => ==(sim), b)
            # band!(ax, bsim.nbon, bsim.low, bsim.upp; alpha=0.3, label=s, color=cols[s])
            lines!(ax, bsim.nbon, bsim.med; label=s, color=cols[s])
        end
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Legend(ga[2,1], ax, orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, focal_sp_range)
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_array_all_meds.png"), fig)

begin
    res = monitored_samplers
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:100:500
    )
    ax1 = Axis(
        gb[1, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax2 = Axis(
        gb[2, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax3 = Axis(
        gb[3, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(ax, b.nbon, b.med .- b.low, b.med .+ b.upp; alpha=0.4, label=s, color=cols[s])
        lines!(ax, b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Legend(ga[2,1], ax, orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, focal_sp_range)
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_array_samplers_med_mean.png"), fig)

begin
    res = monitored_samplers
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:100:500
    )
    ax1 = Axis(
        gb[1, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax2 = Axis(
        gb[2, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax3 = Axis(
        gb[3, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(
            ax,
            b.nbon,
            b.med .- b.lowmed,
            b.med .+ b.uppmed;
            alpha=0.4,
            label=s,
            color=cols[s],
        )
        lines!(ax, b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Legend(ga[2,1], ax, orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, focal_sp_range)
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_array_samplers_med_med.png"), fig)

## Optimized sampling

# Generate BON examples
begin
    Random.seed!(33)
    bons["Focal species range"] = bons["UncertaintySampling"]
    bons["Species richness"] = BON.sample(BON.UncertaintySampling(100), richness_spp)
    bons["Realized interactions"] = BON.sample(
        BON.UncertaintySampling(100), degree_realized
    )
end

# Collect layers
begin
    layers = Dict()
    layers["Focal species range"] = focal_sp_range
    layers["Species richness"] = richness_spp
    layers["Realized interactions"] = degree_realized
end

# Reorder elements for display
_order = Dict(
    "Realized interactions" => 1, "Focal species range" => 2, "Species richness" => 3
)
sort!(monitored_optimized, order(:sampler; by=x -> _order[x]))

# Plot
begin
    res = filter(:sampler => ==("Realized interactions"), sims_optimized)
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:100:500
    )
    ax1 = Axis(
        gb[1, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax2 = Axis(
        gb[2, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax3 = Axis(
        gb[3, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        for sim in unique(res.sim)
            bsim = filter(:sim => ==(sim), b)
            # band!(ax, bsim.nbon, bsim.low, bsim.upp; alpha=0.3, label=s, color=cols[s])
            lines!(ax, bsim.nbon, bsim.med; label=s, color=cols[s])
        end
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, layers[s])
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(
        ga[1, :, Top()], "Optimization layer efficiency"; padding=(0, 0, 5, 0), font=:bold
    )
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_array_optimized_one.png"), fig)

begin
    res = monitored_optimized
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:100:500
    )
    ax1 = Axis(
        gb[1, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax2 = Axis(
        gb[2, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax3 = Axis(
        gb[3, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(ax, b.nbon, b.med - b.low, b.med + b.upp; alpha=0.4, label=s, color=cols[s])
        lines!(ax, b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, layers[s])
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(
        ga[1, :, Top()], "Optimization layer efficiency"; padding=(0, 0, 5, 0), font=:bold
    )
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_array_optimized_med_mean.png"), fig)

begin
    res = monitored_optimized
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:100:500
    )
    ax1 = Axis(
        gb[1, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax2 = Axis(
        gb[2, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    ax3 = Axis(
        gb[3, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(
            ax,
            b.nbon,
            b.med - b.lowmed,
            b.med + b.uppmed;
            alpha=0.4,
            label=s,
            color=cols[s],
        )
        lines!(ax, b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, layers[s])
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(
        ga[1, :, Top()], "Optimization layer efficiency"; padding=(0, 0, 5, 0), font=:bold
    )
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_array_optimized_med_med.png"), fig)

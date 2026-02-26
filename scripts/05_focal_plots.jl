using DrWatson
@quickactivate :NetworkMonitoring

update_theme!(; CairoMakie=(; px_per_unit=2.0))

## Load focal species results

# Use job id to vary parameters
id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
idp = lpad(id, 2, "0")

# Set directory to import results
if !(@isdefined OUTDIR)
    const OUTDIR = "dev" # focal_array or efficiency
end

# Load all results
monitored_types = CSV.read(datadir(OUTDIR, "monitored_types-$idp.csv"), DataFrame)
monitored_types2 = CSV.read(datadir(OUTDIR, "monitored_types2-$idp.csv"), DataFrame)
monitored_spp_all = CSV.read(datadir(OUTDIR, "monitored_spp-$idp.csv"), DataFrame)
monitored_samplers_all = CSV.read(datadir(OUTDIR, "monitored_samplers-$idp.csv"), DataFrame)
monitored_optimized_all = CSV.read(
    datadir(OUTDIR, "monitored_optimized-$idp.csv"), DataFrame
)

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
        @rtransform(:low = :low / :deg, :med = :med / :deg, :upp = :upp / :deg)
        @select(:sim = idp, All())
    end
    return monitored
end
monitored_spp = summarize_monitored(monitored_spp_all)
monitored_samplers = summarize_monitored(monitored_samplers_all)
monitored_optimized = summarize_monitored(monitored_optimized_all)

# Export summarized layers (only for first simulation as example)
if id == 1
    CSV.write(datadir("monitored_spp.csv"), monitored_spp)
    CSV.write(datadir("monitored_samplers.csv"), monitored_samplers)
    CSV.write(datadir("monitored_optimized.csv"), monitored_optimized)
end

# Get species with highest degree
sp = monitored_types.sp[1]
deg = maximum(monitored_types.deg)

# Load layers used for optimization
focal_sp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"))
focal_sp_mask = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_mask-$idp.tiff"))
richness_spp = SDT.SDMLayer(datadir(OUTDIR, "layer_richness_spp-$idp.tiff"))
degree_realized = SDT.SDMLayer(datadir(OUTDIR, "layer_degree_realized-$idp.tiff"))
probsp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_probsp_range-$idp.tiff"))

## Define labels and colors for all plots

# Rename samplers
monitored_samplers.sampler =
    replace.(
        monitored_samplers.sampler,
        "UncertaintySampling" => "Uncertainty Sampling",
        "WeightedBalancedAcceptance" => "Weighted Balanced Acceptance",
        "BalancedAcceptance" => "Balanced Acceptance",
        "SimpleRandomMask" => "Simple Random Mask",
        "SimpleRandom" => "Simple Random",
    )

# Colors
cols = Dict{Any,Any}(
    # Interaction types
    "possible" => Makie.wong_colors()[2],
    "realized" => Makie.wong_colors()[3],
    "detected" => Makie.wong_colors()[4],
    # Samplers
    "Uncertainty Sampling" => Makie.wong_colors()[2],
    "Weighted Balanced Acceptance" => Makie.wong_colors()[3],
    "Simple Random" => Makie.wong_colors()[1],
    "Balanced Acceptance" => :grey,
    "Simple Random Mask" => Makie.wong_colors()[1],
    # Layers
    "Focal species range" => Makie.wong_colors()[2],
    "Species richness" => Makie.wong_colors()[4],
    "Realized interactions" => Makie.wong_colors()[5],
    "Probabilistic range" => Makie.wong_colors()[6],
)
for (sp, col) in zip(unique(monitored_spp.sp), [Makie.wong_colors()[[2, 6, 7]]..., :black])
    cols[sp] = col
end
cols

## Monitored types

begin
    res = monitored_types
    vars = unique(res.var)
    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:25:100
    )
    for v in vars
        b = filter(:var => ==(v), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(b.nbon, b.med; label=v, color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rc, merge=true)
    fig
end
save(plotsdir("focal_types.png"), fig)

# Visualize result
begin
    res = monitored_types2
    vars = unique(res.var)
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:2000:10_000,
    )
    for v in vars
        b = filter(:var => ==(v), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(b.nbon, b.med; label=v, color=cols[v])
        hlines!(unique(b.deg); linestyle=:dash, color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb, merge=true)
    fig
end
save(plotsdir("focal_types2.png"), fig)

## Repeat with 4 species with different degrees

# Visualize result
begin
    res = monitored_spp
    spp = unique(res.sp)
    fig = Figure(; size=(600, 600))
    ax1 = Axis(fig[1, 1]; ylabel="Monitored interactions", xticks=0:100:500, yticks=0:10:60)
    ax2 = Axis(
        fig[2, 1];
        ylabel="Proportion monitored",
        xlabel="Number of sites in BON",
        xticks=0:100:500,
        limits=((nothing, nothing), (0.0, 1.0)),
        yticks=0:0.2:1.0,
    )
    for (i, sp) in enumerate(spp)
        b = filter(:sp => ==(sp), res)
        d = unique(b.deg)[1]
        l = "sp$i: $d int"
        # Monitored int
        band!(ax1, b.nbon, b.low .* b.deg, b.upp .* b.deg; alpha=0.4, label=l)
        lines!(ax1, b.nbon, b.med .* b.deg; label=l)
        hlines!(ax1, d; linestyle=:dash, alpha=0.5)
        # Proportion
        band!(ax2, b.nbon, b.low, b.upp; alpha=0.4, label=l)
        lines!(ax2, b.nbon, b.med; label=l)
    end
    fig[:, end + 1] = Legend(fig, ax1, "Species"; framevisible=false, merge=true)
    fig
end
save(plotsdir("focal_spp.png"), fig)

## Explore variations with different sampler

# Generate BON examples
begin
    Random.seed!(42)
    bons = Dict()
    bons["Uncertainty Sampling"] = BON.sample(BON.UncertaintySampling(100), focal_sp_range)
    bons["Weighted Balanced Acceptance"] = BON.sample(
        BON.WeightedBalancedAcceptance(100), focal_sp_range
    )
    bons["Simple Random"] = BON.sample(BON.SimpleRandom(100), focal_sp_range)
    bons["Balanced Acceptance"] = BON.sample(BON.BalancedAcceptance(100), focal_sp_mask)
    bons["Simple Random Mask"] = BON.sample(BON.SimpleRandom(100), focal_sp_mask)
end

# Plot
begin
    set = ["Uncertainty Sampling", "Balanced Acceptance", "Simple Random Mask"]
    res = @rsubset(monitored_samplers, :sampler in set)
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
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=s, color=cols[s])
        lines!(ax, b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Legend(ga[2,1], ax, orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3], unique(res.sampler))
        heatmap!(a, ifelse(s == "Balanced Acceptance", focal_sp_mask, focal_sp_range))
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    figA = fig
end
save(plotsdir("focal_samplers_mask.png"), fig)

## Optimized sampling

# Generate BON examples
begin
    Random.seed!(33)
    bons["Focal species range"] = bons["Uncertainty Sampling"]
    bons["Species richness"] = BON.sample(BON.UncertaintySampling(100), richness_spp)
    bons["Realized interactions"] = BON.sample(
        BON.UncertaintySampling(100), degree_realized
    )
    bons["Probabilistic range"] = BON.sample(BON.UncertaintySampling(100), probsp_range)
end

# Collect layers
begin
    layers = Dict()
    layers["Focal species range"] = focal_sp_range
    layers["Species richness"] = richness_spp
    layers["Realized interactions"] = degree_realized
    layers["Probabilistic range"] = probsp_range
end

# Reorder elements for display
_order = Dict(
    "Realized interactions" => 1,
    "Focal species range" => 2,
    "Probabilistic range" => 3,
    "Species richness" => 4,
)
sort!(monitored_optimized, order(:sampler; by=x -> _order[x]))

# Plot
begin
    set = [
        "Realized interactions",
        "Focal species range",
        "Species richness",
        "Probabilistic range",
    ]
    res = @rsubset(monitored_optimized, :sampler in set)
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
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=s, color=cols[s])
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
    figB = fig
end
save(plotsdir("focal_optimized.png"), fig)

# Join
begin
    fig = Figure(; size=(650, 1100))
    g1 = GridLayout(fig[1, :])
    g2 = GridLayout(fig[2, :])

    # Figure 1
    res = monitored_samplers
    # Create layouts
    ga = GridLayout(g1[:, 1:3])
    gb = GridLayout(g1[:, end + 1])
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
    ax4 = Axis(
        gb[4, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    hidedecorations!(ax4; label=false)
    # hidespines!(ax4)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=s, color=cols[s])
        lines!(ax, b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="Metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3, ax4], unique(res.sampler))
        heatmap!(a, ifelse(s == "Balanced Acceptance", focal_sp_mask, focal_sp_range))
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Legend
    Legend(
        ga[end + 1, :],
        ax;
        merge=true,
        tellwidth=false,
        tellheight=true,
        nbanks=2,
        framevisible=false,
        labelsize=12.0,
    )

    # Figure 2
    set = [
        "Realized interactions",
        "Focal species range",
        "Probabilistic range",
        "Species richness",
    ]
    res = @rsubset(monitored_optimized, :sampler in set)
    # fig = Figure()
    # Create layouts
    ga = GridLayout(g2[:, 1:3])
    gb = GridLayout(g2[:, end + 1])
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
    ax4 = Axis(
        gb[4, 1]; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10
    )
    # Remove decorations for heatmaps
    hidedecorations!(ax1; label=false)
    hidedecorations!(ax2; label=false)
    hidedecorations!(ax3; label=false)
    hidedecorations!(ax4; label=false)
    # Sampling results
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=s, color=cols[s])
        lines!(ax, b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="Metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    # Heatmaps & BON example
    for (a, s) in zip([ax1, ax2, ax3, ax4], unique(res.sampler))
        heatmap!(a, layers[s])
        scatter!(a, coordinates(bons[s]); markersize=5, color=cols[s], strokewidth=0.5)
        a.ylabel = s
    end
    # Subpanel labels
    Label(
        ga[1, :, Top()], "Optimization layer efficiency"; padding=(0, 0, 5, 0), font=:bold
    )
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Legend
    Legend(
        ga[end + 1, :],
        ax;
        merge=true,
        tellwidth=false,
        tellheight=true,
        nbanks=3,
        framevisible=false,
        labelsize=12.0,
    )

    # Additional labels
    Label(g1[1, :, TopLeft()], "A)"; padding=(0, 0, 5, 0), font=:bold)
    Label(g2[1, :, TopLeft()], "B)"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_joined.png"), fig)
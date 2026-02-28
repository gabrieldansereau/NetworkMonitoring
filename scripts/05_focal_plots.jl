using DrWatson
@quickactivate :NetworkMonitoring

update_theme!(; CairoMakie=(; px_per_unit=2.0))

## Load focal species results

# Use job id to vary parameters
id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
idp = lpad(id, 2, "0")

# Set directory to import results
if !(@isdefined OUTDIR)
    const OUTDIR = "focal_array" # dev(local), focal_array or efficiency
end

# Load test results
monitored_test_all = CSV.read(datadir("monitored_test.csv"), DataFrame)

# Summmarize results not combined previously
function summarize_monitored(df)
    if "layer" in names(df)
        cols = [:set, :sp, :type, :sampler, :layer, :nbon]
    else
        cols = [:set, :sp, :type, :sampler, :nbon]
    end
    monitored = @chain df begin
        groupby(cols)
        @combine(
            :low = quantile(:monitored, 0.05),
            :med = median(:monitored),
            :upp = quantile(:monitored, 0.95),
            :deg = maximum(:deg)
        )
        @rtransform(:low = :low / :deg, :med = :med / :deg, :upp = :upp / :deg)
        @select(:sim = idp, All())
    end
    monitored.sampler =
        replace.(
            monitored.sampler,
            "UncertaintySampling" => "Uncertainty Sampling",
            "WeightedBalancedAcceptance" => "Weighted Balanced Acceptance",
            "BalancedAcceptance" => "Balanced Acceptance",
            "SimpleRandomMask" => "Simple Random Mask",
            "SimpleRandom" => "Simple Random",
        )
    return monitored
end
monitored_test = summarize_monitored(monitored_test_all)

# Define colors for all plots
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
cols
## Monitored types

# Load & summarize results
monitored_types_all = CSV.read(datadir(OUTDIR, "monitored_types-$idp.csv"), DataFrame)
monitored_types2_all = CSV.read(datadir(OUTDIR, "monitored_types2-$idp.csv"), DataFrame)
monitored_types = summarize_monitored(monitored_types_all)
monitored_types2 = summarize_monitored(monitored_types2_all)

# Visualize
fig_types = let
    res = monitored_types
    var = :type
    vals = unique(res[:, var])
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:25:100,
    )
    for v in vals
        b = filter(var => ==(v), res)
        deg = maximum(b.deg)
        band!(b.nbon, b.low * deg, b.upp * deg; alpha=0.4, label=v, color=cols[v])
        lines!(b.nbon, b.med * deg; label=v, color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rc, merge=true)
    fig
end
save(plotsdir("focal_types.png"), fig_types)

# Visualize result
fig_types2 = let
    res = monitored_types2
    var = :type
    vals = unique(res[:, var])
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:2000:10_000,
    )
    for v in vals
        b = filter(var => ==(v), res)
        deg = maximum(b.deg)
        band!(b.nbon, b.low * deg, b.upp * deg; alpha=0.4, label=v, color=cols[v])
        lines!(b.nbon, b.med * deg; label=v, color=cols[v])
        hlines!(deg; linestyle=:dash, color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb, merge=true)
    fig
end
save(plotsdir("focal_types2.png"), fig_types2)

## Repeat with 4 species with different degrees

# Load & summarize results
monitored_spp_all = CSV.read(datadir(OUTDIR, "monitored_spp-$idp.csv"), DataFrame)
monitored_spp = summarize_monitored(monitored_spp_all)
if id == 1
    CSV.write(datadir("monitored_spp.csv"), monitored_spp)
end

# Visualize result
fig_spp = let
    res = monitored_spp
    var = :sp
    vals = unique(res[:, var])
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
    for (i, v) in enumerate(vals)
        b = filter(var => ==(v), res)
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
save(plotsdir("focal_spp.png"), fig_spp)

## Explore variations with different sampler

# Load & summarize results
monitored_samplers_all = CSV.read(datadir(OUTDIR, "monitored_samplers-$idp.csv"), DataFrame)
monitored_samplers = summarize_monitored(monitored_samplers_all)
if id == 1
    CSV.write(datadir("monitored_samplers.csv"), monitored_samplers)
end

# Load layers used for optimization
focal_sp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"))
focal_sp_mask = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_mask-$idp.tiff"))

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
    bons
end

# Plot
fig_samplers = let
    set = ["Uncertainty Sampling", "Weighted Balanced Acceptance", "Simple Random"]
    var = :sampler
    res = filter(var => in(set), monitored_samplers)
    vals = unique(res[:, var])
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:100:500,
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
    for v in vals
        b = filter(var => ==(v), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(ax, b.nbon, b.med; label=v, color=cols[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3], vals)
        heatmap!(a, ifelse(v == "Balanced Acceptance", focal_sp_mask, focal_sp_range))
        scatter!(a, coordinates(bons[v]); markersize=5, color=cols[v], strokewidth=0.5)
        a.ylabel = v
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    figA = fig
end
save(plotsdir("focal_samplers.png"), fig_samplers)

# Plot
fig_mask = let
    set = ["Uncertainty Sampling", "Balanced Acceptance", "Simple Random Mask"]
    var = :sampler
    res = filter(var => in(set), monitored_samplers)
    vals = unique(res[:, var])
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:100:500,
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
    for v in vals
        b = filter(var => ==(v), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(ax, b.nbon, b.med; label=v, color=cols[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3], vals)
        heatmap!(a, ifelse(v == "Uncertainty Sampling", focal_sp_range, focal_sp_mask))
        scatter!(a, coordinates(bons[v]); markersize=5, color=cols[v], strokewidth=0.5)
        a.ylabel = v
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    figA = fig
end
save(plotsdir("focal_samplers_mask.png"), fig_mask)

## Optimized sampling

# Load & summarize results
monitored_optimized_all = CSV.read(
    datadir(OUTDIR, "monitored_optimized-$idp.csv"), DataFrame
)
monitored_optimized = summarize_monitored(monitored_optimized_all)
if id == 1
    CSV.write(datadir("monitored_optimized.csv"), monitored_optimized)
end

# Load layers used for optimization
focal_sp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"))
focal_sp_mask = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_mask-$idp.tiff"))
richness_spp = SDT.SDMLayer(datadir(OUTDIR, "layer_richness_spp-$idp.tiff"))
degree_realized = SDT.SDMLayer(datadir(OUTDIR, "layer_degree_realized-$idp.tiff"))
probsp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_probsp_range-$idp.tiff"))

# Generate BON examples
begin
    Random.seed!(33)
    if !(@isdefined bons)
        bons = Dict()
    end
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
sort!(monitored_optimized, order(:layer; by=x -> _order[x]))

# Plot
fig_optimized = let
    set = [
        "Realized interactions",
        "Focal species range",
        "Species richness",
        "Probabilistic range",
    ]
    var = :layer
    res = filter(var => in(set), monitored_optimized)
    vals = unique(res[:, var])
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:100:500,
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
    for v in vals
        b = filter(var => ==(v), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(ax, b.nbon, b.med; label=v, color=cols[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    example_layers = filter(!=("Probabilistic range"), vals)
    for (a, l) in zip([ax1, ax2, ax3], example_layers)
        heatmap!(a, layers[l])
        scatter!(a, coordinates(bons[l]); markersize=5, color=cols[l], strokewidth=0.5)
        a.ylabel = l
    end
    # Subpanel labels
    Label(
        ga[1, :, Top()],
        "Optimization layer efficiency";
        padding=(0, 0, 5, 0),
        font=:bold,
    )
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    figB = fig
end
save(plotsdir("focal_optimized.png"), fig_optimized)

# Join
fig_joined = let
    ## Define sets to use in both panels
    # Set 1
    res1 = monitored_samplers
    var1 = :sampler
    set1 = [
        "Uncertainty Sampling",
        "Weighted Balanced Acceptance",
        "Simple Random",
        "Balanced Acceptance",
    ]
    # Set 2
    res2 = monitored_optimized
    var2 = :layer
    set2 = [
        "Realized interactions",
        "Focal species range",
        "Probabilistic range",
        "Species richness",
    ]

    # Create main figure elements
    fig = Figure(; size=(650, 1100))
    g1 = GridLayout(fig[1, :])
    g2 = GridLayout(fig[2, :])

    ## Panel 1
    # Define objects
    set = set1
    var = var1
    res = filter(var => in(set), res1)
    vals = unique(res[:, var])

    # Create layouts
    ga = GridLayout(g1[:, 1:3])
    gb = GridLayout(g1[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:100:500,
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
    for v in vals
        b = filter(var => ==(v), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(ax, b.nbon, b.med; label=v, color=cols[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="Metaweb")
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3, ax4], vals)
        heatmap!(a, ifelse(v == "Balanced Acceptance", focal_sp_mask, focal_sp_range))
        scatter!(a, coordinates(bons[v]); markersize=5, color=cols[v], strokewidth=0.5)
        a.ylabel = v
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

    ## Panel 2
    # Define objects
    set = set2
    var = var2
    res = filter(var => in(set), res2)
    vals = unique(res[:, var])

    # Create layouts
    ga = GridLayout(g2[:, 1:3])
    gb = GridLayout(g2[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:100:500,
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
    for v in vals
        b = filter(var => ==(v), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(ax, b.nbon, b.med; label=v, color=cols[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="Metaweb")
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3, ax4], vals)
        heatmap!(a, layers[v])
        scatter!(a, coordinates(bons[v]); markersize=5, color=cols[v], strokewidth=0.5)
        a.ylabel = v
    end

    # Subpanel labels
    Label(
        ga[1, :, Top()],
        "Optimization layer efficiency";
        padding=(0, 0, 5, 0),
        font=:bold,
    )
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

    # Additional labels
    Label(g1[1, :, TopLeft()], "A)"; padding=(0, 0, 5, 0), font=:bold)
    Label(g2[1, :, TopLeft()], "B)"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_joined.png"), fig_joined)

## Range estimation

# Load & summarize results
monitored_estimations_all = CSV.read(
    datadir(OUTDIR, "monitored_estimations-$idp.csv"), DataFrame
)
monitored_estimations = summarize_monitored(monitored_estimations_all)
if id == 1
    CSV.write(datadir("monitored_estimations.csv"), monitored_estimations)
end

# Load layers used for optimization
errors = unique(monitored_estimations.layer)
estimated_ranges = Dict()
for (i, e) in enumerate(errors)
    estimated_ranges[e] = SDT.SDMLayer(
        datadir(OUTDIR, "layer_range_estimations-$idp.tiff"); bandnumber=i
    )
end
estimated_ranges

# Generate BON examples
begin
    Random.seed!(33)
    bons = Dict()
    for e in errors
        bons[e] = BON.sample(BON.BalancedAcceptance(100), estimated_ranges[e])
    end
    bons
end

# Plot
fig_estimation = let
    set = [-0.2, 0.0, 0.2]
    var = :layer
    res = filter(var => in(set), monitored_estimations)
    vals = unique(res[:, var])

    range_over = estimated_ranges[set[1]]
    range_true = estimated_ranges[set[2]]
    range_under = estimated_ranges[set[3]]

    labs = Dict()
    if !(@isdefined cols)
        cols = Dict()
    end
    for (i, s) in enumerate(set)
        cols[s] = Makie.wong_colors()[i]
        labs[s] = string(["Over", "True-", "Under-"][i], s)
    end

    # Create figure
    fig = Figure()
    # Create layouts
    ga = GridLayout(fig[:, 1:3])
    gb = GridLayout(fig[:, end + 1])
    # Create axes
    ax = Axis(
        ga[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:100:500,
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
    for v in vals
        b = filter(var => ==(v), res)
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, color=cols[v], label=labs[v])
        lines!(ax, b.nbon, b.med; label=labs[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    heatmap!(ax1, range_over; colormap=:greys, alpha=0.5)
    heatmap!(ax1, range_true; colormap=:viridis)
    heatmap!(ax2, range_true; colormap=:viridis)
    heatmap!(ax3, range_true; colormap=:viridis, alpha=0.5)
    heatmap!(ax3, range_under; colormap=:viridis)
    for (a, v) in zip([ax1, ax2, ax3], vals)
        scatter!(a, coordinates(bons[v]); markersize=5, strokewidth=0.5, color=cols[v])
        a.ylabel = labs[v]
    end

    # Subpanel labels
    Label(ga[1, :, Top()], "Range estimation"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    fig
end
save(plotsdir("focal_estimation.png"), fig_estimation)

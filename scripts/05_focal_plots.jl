# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module
import BiodiversityObservationNetworks as BON
using BiodiversityObservationNetworks: GI.coordinates

# Use job id to vary parameters
id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
idp = lpad(id, 3, "0")

# Set directory to import results
if !(@isdefined OUTDIR)
    const OUTDIR = "focal_array" # dev(local), focal_array or efficiency
end

# Load & summarize test results
monitored_test_all = CSV.read(datadir("monitored_test.csv"), DataFrame)
monitored_test = summarize_focal(monitored_test_all; id=id)

## Monitored types

# Load & summarize results
monitored_types_all = CSV.read(datadir(OUTDIR, "monitored_types-$idp.csv"), DataFrame)
monitored_types2_all = CSV.read(datadir(OUTDIR, "monitored_types2-$idp.csv"), DataFrame)
monitored_types = summarize_focal(monitored_types_all; id=id)
monitored_types2 = summarize_focal(monitored_types2_all; id=id)

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
        band!(b.nbon, b.low * deg, b.upp * deg; alpha=0.4, label=v, color=colours[v])
        lines!(b.nbon, b.med * deg; label=v, color=colours[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rc, merge=true)
    fig
end
save(plotsdir("supp", "focal_types.png"), fig_types)

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
        band!(b.nbon, b.low * deg, b.upp * deg; alpha=0.4, label=v, color=colours[v])
        lines!(b.nbon, b.med * deg; label=v, color=colours[v])
        hlines!(deg; linestyle=:dash, color=colours[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb, merge=true)
    fig
end
save(plotsdir("supp", "focal_types2.png"), fig_types2)

## Repeat with 4 species with different degrees

# Load & summarize results
monitored_spp_all = CSV.read(datadir(OUTDIR, "monitored_spp-$idp.csv"), DataFrame)
monitored_spp = summarize_focal(monitored_spp_all; id=id)
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
save(plotsdir("supp", "focal_spp.png"), fig_spp)

## Explore variations with different sampler

# Load & summarize results
monitored_samplers_all = CSV.read(datadir(OUTDIR, "monitored_samplers-$idp.csv"), DataFrame)
monitored_samplers = summarize_focal(monitored_samplers_all; id=idp)
if id == 1
    CSV.write(datadir("monitored_samplers.csv"), monitored_samplers)
end

# Load layers used for optimization
focal_sp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"); bandnumber=1)
focal_sp_mask = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"); bandnumber=3)

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
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=colours[v])
        lines!(ax, b.nbon, b.med; label=v, color=colours[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3], vals)
        heatmap!(a, ifelse(v == "Balanced Acceptance", focal_sp_mask, focal_sp_range))
        scatter!(a, coordinates(bons[v]); markersize=5, color=colours[v], strokewidth=0.5)
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
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=colours[v])
        lines!(ax, b.nbon, b.med; label=v, color=colours[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3], vals)
        heatmap!(a, ifelse(v == "Uncertainty Sampling", focal_sp_range, focal_sp_mask))
        scatter!(a, coordinates(bons[v]); markersize=5, color=colours[v], strokewidth=0.5)
        a.ylabel = v
    end
    # Subpanel labels
    Label(ga[1, :, Top()], "Sampler efficiency"; padding=(0, 0, 5, 0), font=:bold)
    Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
    # Show figure
    figA = fig
end
save(plotsdir("supp", "focal_mask.png"), fig_mask)

## Optimized sampling

# Load & summarize results
monitored_optimized_all = CSV.read(
    datadir(OUTDIR, "monitored_optimized-$idp.csv"), DataFrame
)
monitored_optimized = summarize_focal(monitored_optimized_all; id=id)
if id == 1
    CSV.write(datadir("monitored_optimized.csv"), monitored_optimized)
end

# Load layers used for optimization
focal_sp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"); bandnumber=1)
focal_sp_mask = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"); bandnumber=3)
probsp_range = SDT.SDMLayer(datadir(OUTDIR, "layer_sp_range-$idp.tiff"); bandnumber=2)
richness_spp = SDT.SDMLayer(datadir(OUTDIR, "layer_richness-$idp.tiff"); bandnumber=1)
degree_realized = SDT.SDMLayer(datadir(OUTDIR, "layer_degree-$idp.tiff"); bandnumber=1)

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
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=colours[v])
        lines!(ax, b.nbon, b.med; label=v, color=colours[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    # axislegend(ax; position=:lt, merge=true, labelsize=14)
    Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=2)
    # Heatmaps & BON example
    example_layers = filter(!=("Probabilistic range"), vals)
    for (a, l) in zip([ax1, ax2, ax3], example_layers)
        heatmap!(a, layers[l])
        scatter!(a, coordinates(bons[l]); markersize=5, color=colours[l], strokewidth=0.5)
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
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=colours[v])
        lines!(ax, b.nbon, b.med; label=v, color=colours[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="Metaweb")
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3, ax4], vals)
        heatmap!(a, ifelse(v == "Balanced Acceptance", focal_sp_mask, focal_sp_range))
        scatter!(a, coordinates(bons[v]); markersize=5, color=colours[v], strokewidth=0.5)
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
        band!(ax, b.nbon, b.low, b.upp; alpha=0.4, label=v, color=colours[v])
        lines!(ax, b.nbon, b.med; label=v, color=colours[v])
    end
    hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey, label="Metaweb")
    # Heatmaps & BON example
    for (a, v) in zip([ax1, ax2, ax3, ax4], vals)
        heatmap!(a, layers[v])
        scatter!(a, coordinates(bons[v]); markersize=5, color=colours[v], strokewidth=0.5)
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

# Use job id to vary parameters
# id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
# idp = lpad(id, 3, "0")

# Default to efficiency results for illustration
const OUTDIR = "efficiency" # dev(local), focal_array or efficiency

# Load & summarize results
monitored_estimations_all = CSV.read(
    datadir(OUTDIR, "monitored_estimations-$idp.csv"), DataFrame
)
monitored_missings = filter(:monitored => ismissing, monitored_estimations_all)
filter!(:monitored => !ismissing, monitored_estimations_all)
monitored_estimations = summarize_focal(
    monitored_estimations_all; id=id, confint=true, α=0.10
)
@rtransform!(monitored_estimations, :degmax = :degmax / :deg)
rename!(monitored_estimations, :degmax => :pmax)
if id == 1
    CSV.write(datadir("monitored_estimations.csv"), monitored_estimations)
end

# Load layers used for optimization
errors = string.(unique(monitored_estimations.layer))
offsets = unique(monitored_estimations.offset)
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
    bons_adj = Dict()
    for (i, e) in enumerate(errors)
        n = 100
        n_adj = round(Int, ((1 + offsets[i]) * n))
        bons[e] = BON.sample(BON.BalancedAcceptance(n), estimated_ranges[e])
        bons_adj[e] = BON.sample(BON.BalancedAcceptance(n_adj), estimated_ranges[e])
    end
    bons
end

# Add effort adjustment
nbon_ref = maximum(@rsubset(monitored_estimations, :offset == 0.0).nbon)
nbon_max = nbon_ref
@chain monitored_estimations begin
    @rtransform!(:nmax = round(Int, nbon_ref * (1 + :offset)))
    @rtransform!(:neff = :nbon * nbon_max / :nmax)
end

# Plot
fig_estimation = begin
    function plot_focal(;
        adjust_effort=false, adjust_n=false, option=:integral, n=300, p=0.95, pmax=false
    )
        set = ["Over-0.50", "True-0.00", "Under-0.50"]
        var = :layer
        res = disallowmissing(filter(var => in(set), monitored_estimations))
        vals = unique(res[:, var])

        # Display elements
        show_lines = true
        show_scatter = true
        show_eff = false
        show_sat = true
        show_int = false
        adjust_effort = adjust_effort
        adjust_n = adjust_n

        # Adjust sampling effort
        _bons = bons_adj
        if adjust_effort
            @rsubset!(res, :nbon <= :nmax)
        end
        if adjust_n
            @rtransform!(res, :nbon = :neff)
        end

        # Replace set for illustration when overestimation is not available.
        replaced = false
        if !(all([s in errors for s in set]))
            _errs = filter(startswith("Over-"), errors)
            _max =
                parse.(Float64, replace.(_errs, "Over-" => "")) |>
                maximum |>
                s -> @sprintf("%.2f", s)
            _orig_set = set
            set = ["Over-$_max", "True-0.00", "Under-$_max"]
            if all([s in errors for s in set])
                @warn "Requested set not available in exported layers. Replacing by closest set: $set"
                replaced = true
                res = disallowmissing(filter(var => in(set), monitored_estimations))
                vals = unique(res[:, var])
            else
                @error "Requested set $set not available in exported layers"
            end
        end

        # Set layers
        range_over = estimated_ranges[set[1]]
        range_true = estimated_ranges[set[2]]
        range_under = estimated_ranges[set[3]]

        # Set colours
        if !(@isdefined colours)
            colours = Dict()
        end
        for (i, s) in enumerate(set)
            colours[s] = Makie.wong_colors()[i]
        end

        # Create figure
        fig = Figure(; size=(600, 500))
        # Create layouts
        ga = GridLayout(fig[:, 1:3])
        gb = GridLayout(fig[:, end + 1])
        # Create axes
        xlim = 500
        ax = Axis(
            ga[1, 1:3];
            xlabel="Sites in BON",
            ylabel="Monitored interactions",
            xticks=0:100:500,
            limits=((nothing, 1.02 * xlim), (nothing, nothing)),
        )
        ax0 = Axis(ga[2, 2:3]; xlabel="Efficiency")
        yopts = (; aspect=1, yaxisposition=:right, ylabelrotation=1.5pi, ylabelsize=10)
        ax1 = Axis(gb[1, 1]; yopts...)
        ax2 = Axis(gb[2, 1]; yopts...)
        ax3 = Axis(gb[3, 1]; yopts...)
        # Highlight replaced set values
        if replaced
            ax1.ylabelcolor = :red
            ax3.ylabelcolor = :red
        end
        # Remove decorations for heatmaps
        hidedecorations!(ax1; label=false)
        hidedecorations!(ax2; label=false)
        hidedecorations!(ax3; label=false)

        # Sampling results
        for (i, v) in enumerate(vals)
            b = filter(var => ==(v), res)
            # Get the saturation parameter for the curve
            pm = pmax ? unique(b.pmax)[1] : 1.0
            eff_a = efficiency_gridsearch(b.nbon, b.med, pm; f=exp)
            @info "$v: a = $eff_a"
            eff_a_low = efficiency_gridsearch(b.nbon, b.confint_low, pm; f=exp)
            eff_a_upp = efficiency_gridsearch(b.nbon, b.confint_upp, pm; f=exp)
            # Get the efficiencies for comparison
            nv = n isa Dict ? n[v] : n
            eff = efficiency(b.nbon, b.med; f=exp, pmax=pm, option=option, n=nv, p=p)
            eff_low = efficiency(
                b.nbon, b.confint_low; f=exp, pmax=pm, option=option, n=nv, p=p
            )
            eff_upp = efficiency(
                b.nbon, b.confint_upp; f=exp, pmax=pm, option=option, n=nv, p=p
            )
            if v == "True-0.00"
                global _nbon = b.nbon
                global _med = b.med
            end
            # Display results
            lab = ifelse(show_eff, "$v, eff=$(@sprintf("%.0f", eff))", v)
            band!(ax, b.nbon, b.low, b.upp; alpha=0.4, color=colours[v], label=lab)
            if show_sat
                lines!(
                    ax,
                    1:xlim,
                    saturation(eff_a, pm)(1:xlim);
                    color=colours[v],
                    linestyle=:dash,
                    linewidth=1.5,
                    alpha=0.7,
                )
            end
            if show_int
                lines!(
                    ax,
                    1:xlim,
                    saturation(eff_a_low, pm)(1:xlim);
                    color=colours[v],
                    linestyle=:dot,
                    linewidth=1.5,
                )
                lines!(
                    ax,
                    1:xlim,
                    saturation(eff_a_upp, pm)(1:xlim);
                    color=colours[v],
                    linestyle=:dot,
                    linewidth=1.5,
                )
            end
            if show_lines
                lines!(ax, b.nbon, b.med; label=lab)
            end
            if show_scatter
                scatter!(ax, b.nbon, b.med; label=lab, markersize=5)
            end
            rangebars!(ax0, [i], [eff_low], [eff_upp]; direction=:x, color=colours[v])
            if v == "True-0.00"
                vlines!(ax0, [eff]; color=colours[v], linestyle=:dash, alpha=0.5)
                band!(
                    ax0,
                    [eff_low, eff_upp],
                    [0.0, 0.0],
                    [4.0, 4.0];
                    color=colours[v],
                    alpha=0.4,
                )
            end
            scatter!(ax0, [eff], [i]; color=colours[v])
        end
        hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey)
        Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=3)
        # Heatmaps & BON example
        heatmap!(ax1, range_over; colormap=:greys, alpha=0.5)
        heatmap!(ax1, range_true; colormap=:viridis)
        heatmap!(ax2, range_true; colormap=:viridis)
        heatmap!(ax3, range_true; colormap=:viridis, alpha=0.5)
        heatmap!(ax3, range_under; colormap=:viridis)
        for (a, v) in zip([ax1, ax2, ax3], vals)
            scatter!(a, coordinates(_bons[v]); markersize=5, strokewidth=0.5, color=colours[v])
            a.ylabel = v
        end

        # Fix efficiency panel
        limits!(ax0, (nothing, nothing), (0.5, 3.5))
        hideydecorations!(ax0)
        ax0.yreversed = true
        # Adapt panel for options
        adj_n = adjust_n ? "(n ajusted)" : ""
        adj_e = adjust_effort ? "(effort-adjusted)" : ""
        pm_t = pmax ? "(pmax = true)" : ""
        if option == :integral
            ax0.xlabel = "Efficiency integral $(pm_t)$(adj_n)$(adj_e)"
        elseif option == :n_at_p
            ax0.xlabel = "Number of sites at p = $p $(pm_t)$(adj_n)$(adj_e)"
        elseif option == :p_at_n
            if n isa Dict
                ax0.xlabel = "Proportion at maximum n"
            else
                ax0.xlabel = "Proportion at n = $n $(pm_t)$(adj_n)$(adj_e)"
            end
        elseif option == :integral_at_n
            ax0.xlabel = "Efficiency integral at n = $n $(pm_t)$(adj_n)$(adj_e)"
        elseif option == :a
            ax0.xlabel = "Parameter a"
        end

        # Subpanel labels
        Label(
            ga[1, :, Top()],
            "Range estimation, id=$idp";
            padding=(0, 0, 5, 0),
            font=:bold,
        )
        Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
        # Show figure
        save(plotsdir("focal_ranges.png"), fig)
        return fig
    end
    plot_focal(; adjust_effort=false, adjust_n=false)
end

## Adjust n

# Adjusting n at 300
plot_focal(; adjust_effort=false, adjust_n=false)
plot_focal(; adjust_effort=false, adjust_n=true)

# Integral at n
plot_focal(; adjust_effort=false, adjust_n=false, option=:integral)
plot_focal(; adjust_effort=false, adjust_n=false, option=:integral_at_n, n=500)
plot_focal(; adjust_effort=false, adjust_n=true, option=:integral_at_n, n=500)

# Evaluating n at p = 0.95
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_p, p=0.95)
plot_focal(; adjust_effort=false, adjust_n=true, option=:n_at_p, p=0.95)

# Evaluating p at n
plot_focal(; adjust_effort=false, adjust_n=false, option=:p_at_n, n=300)
plot_focal(;
    adjust_effort=false,
    adjust_n=false,
    option=:p_at_n,
    n=Dict("Over-0.50" => 450, "True-0.00" => 300, "Under-0.50" => 150),
)

## pmax

# integral
plot_focal(; adjust_effort=false, adjust_n=false, pmax=false)
plot_focal(; adjust_effort=false, adjust_n=false, pmax=true)

# integral at n
plot_focal(; adjust_effort=false, adjust_n=false, option=:integral_at_n, n=5000, pmax=false)
plot_focal(; adjust_effort=false, adjust_n=false, option=:integral_at_n, n=5000, pmax=true)
plot_focal(; adjust_effort=false, adjust_n=false, option=:integral_at_n, n=500, pmax=true)
plot_focal(; adjust_effort=false, adjust_n=false, option=:integral_at_n, n=300, pmax=true)

# p at n
plot_focal(; adjust_effort=false, adjust_n=false, option=:p_at_n, n=300, pmax=false)
plot_focal(; adjust_effort=false, adjust_n=false, option=:p_at_n, n=300, pmax=true)
plot_focal(; adjust_effort=false, adjust_n=false, option=:p_at_n, n=500, pmax=true)

# n at p
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_p, p=0.80, pmax=false)
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_p, p=0.80, pmax=true)
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_p, p=0.50, pmax=true)

# a - don't use a with pmax=true, only an option for convenience
plot_focal(; adjust_effort=false, adjust_n=false, option=:a, pmax=false)
plot_focal(; adjust_effort=false, adjust_n=false, option=:a, pmax=true)

p = 0.8
pmax = 0.8947368421052632
a_true = 333.3126399012809
a_over = 497.2639203889265
a_under = 356.7687128894338
a_under_pmax = 308.2999298366703

a_under / a_under_pmax
a_under_pmax / a_under
a_under_pmax / p

efficiency_n_at_p(a_true, p, 1.0)
efficiency_n_at_p(a_over, p, 1.0)
# références
efficiency_n_at_p(a_under, p, 1.0)
# plus haut, juste mais pas bien étalonné
efficiency_n_at_p(a_under_pmax, p * pmax, pmax)
# bien étalonné, mais trop bas
efficiency_n_at_p(a_under_pmax, p, 1.0)
# revient au même
efficiency_n_at_p(a_under_pmax / pmax, p * pmax, pmax)
# bien étalonné, juste que ce soit plus bas
# est-ce que c'est valide dans tous les cas?
efficiency_n_at_p(a_under_pmax, p * pmax, pmax) / pmax
# équivalent
efficiency_n_at_p(a_under_pmax, p * pmax, pmax) / p
# effet trop important? justifiable? quand même moins que over
a_under / a_under_pmax
a_under / a_under_pmax * efficiency_n_at_p(a_under_pmax, p * pmax, pmax)
# reverse-engineer de valeur sans correction! Genre de produit croisé? Justifiable ?

# autre alternative avec

# n at pmax
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_pmax, p=0.8, pmax=false)
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_pmax, p=0.8, pmax=true)
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_pmax2, p=0.8, pmax=true)
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_pmax3, p=0.8, pmax=true)
plot_focal(; adjust_effort=false, adjust_n=false, option=:n_at_pmax4, p=0.8, pmax=false)

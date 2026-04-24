# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module

# Load data
effs_samplers = CSV.read(datadir("efficiency_samplers.csv"), DataFrame)
effs_optimized = CSV.read(datadir("efficiency_optimized.csv"), DataFrame)
effs_species = CSV.read(datadir("efficiency_species.csv"), DataFrame)

# Remove Simple Random Mask
@rsubset!(effs_samplers, :variable != "Simple Random Mask")

# Convert outliers - Efficiency should not go beyond 10,000 sites (maximum in landscape)
# when efficiency is measured as the number of sites to reach 80% of interactions
for effs in [effs_samplers, effs_optimized, effs_species]
    lim = 10_000.0
    @rtransform! effs begin
        :eff = :eff > lim ? lim : :eff
        :eff_low = :eff_low > lim ? lim : :eff_low
        :eff_upp = :eff > lim ? lim : :eff_upp
    end
end

# Combine for convenience
effs_combined = vcat(effs_samplers, effs_optimized; cols=:union)

# Define sorting order across figures
sortedsamplers = [
    "Uncertainty Sampling",
    "Weighted Balanced Acceptance",
    "Simple Random",
    "Balanced Acceptance",
    "Simple Random Mask",
]
sortedlayers = [
    "Focal species range",
    "Realized interactions",
    "Probabilistic range",
    "Species richness",
]
sortedlayout = [sortedsamplers..., sortedlayers...]
sorteddict = Dict(v => i for (i, v) in enumerate(sortedlayout))

# Sort dataframes
sort!(effs_samplers, order(:variable; by=x -> sorteddict[x]))
sort!(effs_optimized, order(:variable; by=x -> sorteddict[x]))
sort!(effs_combined, order(:variable; by=x -> sorteddict[x]))

## Within-simulation comparisons

# Separate results per simulation
compsdict = Dict(
    "Uncertainty Sampling" => "US",
    "Weighted Balanced Acceptance" => "WBA",
    "Simple Random" => "RS",
    "Balanced Acceptance" => "BA",
    "Simple Random Mask" => "SRM",
    "Focal species range" => "FR",
    "Realized interactions" => "RI",
    "Probabilistic range" => "PR",
    "Species richness" => "SR",
)
set = unique(effs_combined.variable)
within_comps_all = comparewithin(
    effs_combined, set; labels=compsdict, f=(n1, n2) -> n1 / n2
)

# Let's flip a comparison for illustration
all_comps = unique(within_comps_all.variable)
toflip = all_comps
flipthatcomp!(within_comps_all, toflip; f=n -> 1 / n)

# Count positive and negative comparisons per set and variable (across simulations/replicates)
unique_comps_all = @chain within_comps_all begin
    @groupby :set :variable
    @combine(
        :n = length(:value),
        :sign_pos = count(>(1), :value) / length(:value),
        :sign_neg = count(<=(0), :value) / length(:value),
        :positive = count(==("positive"), :overlap) / length(:overlap),
        :negative = count(==("negative"), :overlap) / length(:overlap),
        :overlap = count(==("overlap"), :overlap) / length(:overlap),
    )
    @transform :prop_sim = :n ./ maximum(:n)
    stack(
        [:sign_pos, :sign_neg, :positive, :negative, :overlap];
        variable_name=:countmeasure,
        value_name=:count,
    )
    @rtransform :label = "$(round(Int, :count *100)) %"
    @rtransform :label = (:count > 0.0 && :label == "0 %") ? "< 1 %" : :label
end

# Select a reduced number of comparisons
reducedcomps_dict = Dict(
    # Samplers
    "ΔRS_US" => "Simple Random",
    "ΔWBA_US" => "Weighted Balanced Acceptance",
    "ΔBA_US" => "Balanced Acceptance",
    "ΔSRM_US" => "Simple Random Mask",
    # Layers
    "ΔRI_FR" => "Realized Interactions",
    "ΔSR_FR" => "Species Richness",
    "ΔPR_FR" => "Probabilistic Range",
)
# Select them
within_comps = @rsubset within_comps_all :variable in collect(keys(reducedcomps_dict))
unique_comps = @rsubset res_summary_all :variable in collect(keys(reducedcomps_dict))
# Rename simply
@rtransform! within_comps :variable = reducedcomps_dict[:variable]
@rtransform! unique_comps :variable = reducedcomps_dict[:variable]

## Plot the comparison

# Visualize
begin
    # Select results for comparison
    res_comps = within_comps
    res_summary = @rsubset(
        unique_comps, :countmeasure in ["positive", "negative", "overlap"]
    )
    sortedcomps = unique(res_summary.variable)

    # Set colour palette
    pal = Dict(
        "negative" => Makie.wong_colors()[2],
        "overlap" => Makie.wong_colors()[4],
        "positive" => Makie.wong_colors()[3],
    )
    scl = scales(; Color=(; palette=[k => v for (k, v) in pal]))
end;
let
    Random.seed!(42)

    # Figure
    f = Figure(; size=(850, 450))
    g1 = GridLayout(f[1:3, 1:4])
    g2 = GridLayout(f[4:6, 1:4])
    g3 = GridLayout(f[1:3, (end + 1):(end + 2)])
    g4 = GridLayout(f[4:end, (end - 1):end])

    # Main panels
    d1 = @rsubset(res_comps, :set == "samplers")
    d2 = @rsubset(res_comps, :set == "layers")
    m = mapping(:variable => sorter(sortedcomps) => "", :value => log2; color=:overlap)
    rains = visual(
        RainClouds;
        markersize=5,
        jitter_width=0.3,
        plot_boxplots=false,
        clouds=nothing,
        orientation=:horizontal,
    )
    vline = mapping([0.0]) * visual(VLines; linestyle=:dash)
    # Tweak axis
    xlog2f = vs -> [rich("2", superscript("$(Int(v))")) for v in vs]
    xticks = [-6, -4, -2, 0, 2, 4]
    xaxis = (; xticks=xticks, xtickformat=xlog2f)
    xlab1 = "Efficiency compared to reference (Uncertainty Sampling)"
    xlab2 = "Efficiency compared to reference (Focal Range)"
    # Draw figures
    fg1 = draw!(g1, data(d1) * m * rains + vline, scl; axis=(; xlabel=xlab1, xaxis...))
    fg2 = draw!(g2, data(d2) * m * rains + vline, scl; axis=(; xlabel=xlab2, xaxis...))
    linkxaxes!(fg1..., fg2...)
    # Add labels
    pad = (-150, 0, 10, 0)
    Label(g1[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(g2[1, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)

    # Summary panels
    d3 = @rsubset(res_summary, :set == "samplers")
    d4 = @rsubset(res_summary, :set == "layers")
    ax3 = Axis(g3[1, 1]; xticks=([0.5, 1.5, 2.5], ["Negative", "Overlap", "Positive"]))
    ax4 = Axis(g4[1, 1]; xticks=([0.5, 1.5, 2.5], ["Negative", "Overlap", "Positive"]))
    m34 = mapping(
        :variable => sorter(sortedcomps),
        [1];
        stack=:countmeasure => sorter(["negative", "overlap", "positive"]),
        color=:countmeasure,
        bar_labels=:label => verbatim,
    )
    v3 = visual(
        BarPlot;
        direction=:x,
        label_position=:center,
        label_color=:white,
        label_font=:bold,
        label_size=14,
        alpha=0.85,
    )
    v4 = visual(
        BarPlot;
        direction=:x,
        label_position=:center,
        label_color=:white,
        label_font=:bold,
        label_size=14,
        alpha=0.85,
    )
    draw!(ax3, data(d3) * m34 * v3, scl)
    draw!(ax4, data(d4) * m34 * v4, scl)

    hideydecorations!(ax3)
    hidexdecorations!(ax3; ticklabels=false, ticks=false)
    hidespines!(ax3)

    hideydecorations!(ax4)
    hidexdecorations!(ax4; ticklabels=false, ticks=false)
    hidespines!(ax4)

    pad = (0, 0, 10, 0)
    Label(g3[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)
    Label(g4[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)

    save(plotsdir("efficiency_comparison.png"), current_figure())
    f
end

## Within-simulation - All comparisons

# Visualize
sortedcomps = [
    # Samplers
    "ΔBA_SRM"
    "ΔUS_BA"
    "ΔUS_SRM"
    "ΔUS_RS"
    "ΔUS_WBA"
    "ΔWBA_RS"
    # Layers
    "ΔFR_PR"
    "ΔFR_SR"
    "ΔPR_SR"
    "ΔRI_FR"
    "ΔRI_PR"
    "ΔRI_SR"
]
let res_comps = within_comps_all, res_summary = res_summary_all
    Random.seed!(42)
    d0 = @rsubset(res_comps, :variable in sortedcomps)
    u0 = @rsubset(res_summary, :variable in sortedcomps)

    # Figure & grid
    f = Figure(; size=(700, 450))
    g1 = GridLayout(f[1:4, 1:3])
    g2 = GridLayout(f[5:10, 1:3])
    g3 = GridLayout(f[1:4, end + 1])
    g4 = GridLayout(f[5:end, end])

    # Main panels
    d1 = @rsubset(d0, :set == "samplers")
    d2 = @rsubset(d0, :set == "layers")
    m = mapping(
        :variable => sorter(sortedcomps) => "comparison",
        :value => "Efficiency difference";
        color=:value => (x -> x <= 0.0),
    )
    rains = visual(
        RainClouds;
        markersize=5,
        jitter_width=0.4,
        plot_boxplots=false,
        clouds=nothing,
        orientation=:horizontal,
    )
    vline = mapping([0.0]) * visual(VLines; linestyle=:dash)
    # xlog2f = vs -> [rich("2", superscript("$(v)")) for v in vs]
    # xlog2 = (; axis=(; xtickformat=xlog2f))
    # xaxis = (; axis=(; xticks=0:2:10))
    fg1 = draw!(g1, data(d1) * m * rains + vline;)
    fg2 = draw!(g2, data(d2) * m * rains + vline;)
    linkxaxes!(fg1..., fg2...)
    pad = (-100, 0, 10, 0)
    Label(g1[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(g2[1, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)

    # Summary panels
    d3 = @rsubset(u0, :set == "samplers")
    d4 = @rsubset(u0, :set == "layers")
    ax3 = Axis(g3[1, 1]; xticks=([0.5, 1.5], ["Negative", "Positive"]))
    ax4 = Axis(g4[1, 1]; xticks=([0.5, 1.5], ["Negative", "Positive"]))
    m34 = mapping(
        :variable => sorter(sortedcomps),
        [1];
        stack=:countmeasure => sorter(["count_neg", "count_pos"]),
        color=:countmeasure => sorter(["count_pos", "count_neg"]),
        bar_labels=:label => verbatim,
    )
    v3 = visual(
        BarPlot;
        direction=:x,
        label_position=:center,
        label_color=:white,
        label_font=:bold,
        label_size=14,
        alpha=0.85,
    )
    v4 = visual(
        BarPlot;
        direction=:x,
        label_position=:center,
        label_color=:white,
        label_font=:bold,
        label_size=14,
        alpha=0.85,
    )
    draw!(ax3, data(d3) * m34 * v3)
    draw!(ax4, data(d4) * m34 * v4)

    hideydecorations!(ax3)
    hidexdecorations!(ax3; ticklabels=false, ticks=false)
    hidespines!(ax3)

    hideydecorations!(ax4)
    hidexdecorations!(ax4; ticklabels=false, ticks=false)
    hidespines!(ax4)

    pad = (0, 0, 10, 0)
    Label(g3[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)
    Label(g4[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)

    save(plotsdir("supp", "efficiency_comparison_all.png"), current_figure())
    f
end
## Distribution figures

# Violin
begin
    Random.seed!(42) # for jitter
    eff = :eff => "efficiency"
    layer =
        mapping(:variable => sorter(reverse(sortedlayout)) => "", eff; color=:variable) *
        visual(
            RainClouds;
            markersize=5,
            jitter_width=0.25,
            plot_boxplots=false,
            orientation=:horizontal,
            # dodge_gap=0.2, # not working
            gap=0.1,
            cloud_width=1.0,
        )
    f1 = data(effs_samplers) * layer
    f2 = data(effs_optimized) * layer
    f = Figure(; size=(700, 700))
    sc = scales(; Color=(; palette=colourpal))
    ylog2f = (; ytickformat=vs -> [rich("2", superscript("$(Int(v))")) for v in vs])
    ylog2 = (; axis=(;))
    fg1 = draw!(f[1, 1], f1, sc; ylog2...)
    fg2 = draw!(f[2, 1], f2, sc; ylog2...)
    linkyaxes!(fg1..., fg2...)
    pad = (-45, 0, 20, 0)
    Label(f[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(f[2, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)
    save(plotsdir("supp", "efficiency_distribution.png"), f)
    f
end

# Efficiency given species degree
begin
    Random.seed!(42)
    ranknames = renamer(4 => "0.07", 3 => "0.33", 2 => "0.66", 1 => "1.0")
    ranks = :rank => ranknames => "Within-simulation Degree Percentile Rank"
    sppalette = from_continuous(cgrad(:viridis; rev=false))
    f = Figure(; size=(650, 650))
    fg1 = draw!(
        f[1, 1],
        data(effs_species) *
        mapping(ranks, eff; color=ranks) *
        visual(
            RainClouds;
            markersize=5,
            jitter_width=0.1,
            plot_boxplots=false,
            # orientation=:horizontal,
        ),
        scales(; Color=(; palette=sppalette));
        axis=(;),
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
    save(plotsdir("supp", "efficiency_distribution_species.png"), f)
    f
end

## Efficiency & Occupancy

# Scatter & smooth with independent subfigures
begin
    occ = :occ => "occupancy"
    layout = mapping(occ, eff; color=:variable) * (visual(Scatter) + linear())
    scl = scales(; Color=(; palette=colourpal))
    legend = (; position=:bottom)
    f1 = data(effs_samplers) * layout * mapping(; col=:variable => sorter(sortedlayout))
    f2 = data(effs_optimized) * layout * mapping(; col=:variable => sorter(sortedlayout))
    f = Figure(; size=(900, 500))
    fg1 = draw!(f[1, 1:4], f1, scl; ylog2...)
    fg2 = draw!(f[2, 1:4], f2, scl; ylog2...)
    linkaxes!(fg1..., fg2...)
    pad = (-50, 0, 30, 0)
    Label(f[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(f[2, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)
    save(plotsdir("supp", "efficiency_occupancy.png"), f)
    f
end

# Scatter & smooth with percentile rank
begin
    sppalette = from_continuous(cgrad(:viridis; rev=false))
    layout = mapping(occ, eff; color=:variable) * (visual(Scatter) + linear())
    f = draw(
        data(effs_species) * layout * mapping(; color=ranks, col=ranks),
        scales(; Color=(; palette=sppalette));
        figure=(; size=(800, 300)),
        legend=legend,
        axis=(;),
    )
    save(plotsdir("supp", "efficiency_occupancy_species.png"), f)
    f
end
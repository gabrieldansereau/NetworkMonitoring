# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module

# Load data
effs_samplers = CSV.read(datadir("efficiency_samplers.csv"), DataFrame)
effs_optimized = CSV.read(datadir("efficiency_optimized.csv"), DataFrame)
effs_species = CSV.read(datadir("efficiency_species.csv"), DataFrame)

# Remove Simple Random Mask
@rsubset!(effs_samplers, :variable != "Simple Random Mask")

# Combine for convenience
effs_combined = vcat(effs_samplers, effs_optimized; cols=:union)

## Efficiency only

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

## Efficiency & Occupancy with species degree

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

## Within-simulation comparison

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
within_combined_dif = comparewithin(
    effs_combined,
    set;
    labels=compsdict,
    # f=(n, n2) -> efficiency_difference(n, n2; k=10_000),
    f=(n, n2) -> n - n2,
)

# Let's flip a comparison for illustration
toflip = ["ΔFR_RI"]
flipthatcomp!(within_combined_dif, toflip)

# Count number of positive and negative comparisons
unique_df = @chain within_combined_dif begin
    @groupby(:set, :variable)
    @combine(
        :count_pos = count(>(0), :value) / length(:value),
        :count_neg = count(<=(0), :value) / length(:value)
    )
    stack([:count_pos, :count_neg]; variable_name=:countmeasure, value_name=:count)
    @rtransform(:label = "$(round(Int, :count *100)) %")
    @rtransform(:label = (:count > 0.0 && :label == "0 %") ? "< 1 %" : :label)
end

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
let d = within_combined_dif, u = unique_df
    Random.seed!(42)
    d0 = @rsubset(d, :variable in sortedcomps)
    u0 = @rsubset(u, :variable in sortedcomps)

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

## Reduced number of comparison

# Reduce number of comparisons
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
within_combined_dif2 = comparewithin(
    effs_combined,
    set;
    to=["Uncertainty Sampling", "Focal species range"],
    labels=compsdict,
    # f=(n, n2) -> efficiency_difference(n, n2; k=10_000),
    f=(n, n2) -> n - n2,
)
flipthatcomp!(within_combined_dif2, unique(within_combined_dif2.variable))
@rtransform!(within_combined_dif2, :variable = reducedcomps_dict[:variable])

# Count number of positive and negative comparisons
unique_df2 = @chain within_combined_dif2 begin
    @groupby(:set, :variable)
    @combine(
        :count_pos = count(>(0), :value) / length(:value),
        :count_neg = count(<=(0), :value) / length(:value)
    )
    stack([:count_pos, :count_neg]; variable_name=:countmeasure, value_name=:count)
    @rtransform(:label = "$(round(Int, :count *100)) %")
    @rtransform(:label = (:count > 0.0 && :label == "0 %") ? "< 1 %" : :label)
end

# Visualize
let d = within_combined_dif2, u = unique_df2, sortedcomps = unique(u.variable)
    Random.seed!(42)

    # Figure
    f = Figure(; size=(750, 450))
    g1 = GridLayout(f[1:3, 1:5])
    g2 = GridLayout(f[4:6, 1:5])
    g3 = GridLayout(f[1:3, (end + 1):(end + 2)])
    g4 = GridLayout(f[4:end, (end - 1):end])

    # Main panels
    d1 = @rsubset(d, :set == "samplers")
    d2 = @rsubset(d, :set == "layers")
    m = mapping(
        :variable => sorter(sortedcomps) => "",
        :value => "Efficiency compared to reference (Uncertainty Sampling)";
        color=:value => (x -> x <= 0.0),
    )
    rains = visual(
        RainClouds;
        markersize=5,
        jitter_width=0.3,
        plot_boxplots=false,
        clouds=nothing,
        orientation=:horizontal,
    )
    vline = mapping([0.0]) * visual(VLines; linestyle=:dash)

    xlog2f = vs -> [rich("2", superscript("$(v)")) for v in vs]
    xlog2 = (; axis=(; xtickformat=xlog2f))
    # xaxis = (; axis=(; xticks=0:2:10))
    fg1 = draw!(g1, data(d1) * m * rains + vline;)
    fg2 = draw!(
        g2,
        data(d2) * m * rains + vline;
        axis=(; xlabel="Efficiency compared to reference (Focal Range)"),
    )
    linkxaxes!(fg1..., fg2...)
    pad = (-150, 0, 10, 0)
    Label(g1[1, 1, Top()], "A) Samplers"; halign=:left, font=:bold, padding=pad)
    Label(g2[1, 1, Top()], "B) Optimization Layers"; halign=:left, font=:bold, padding=pad)

    # Summary panels
    d3 = @rsubset(u, :set == "samplers")
    d4 = @rsubset(u, :set == "layers")
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

    save(plotsdir("efficiency_comparison_reduced.png"), current_figure())
    f
end

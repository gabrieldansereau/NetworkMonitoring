## Extra figures

# Run main script first
include("09_ranges_plots.jl")

## Add sampling effort over area

# Add area and sampling effort
@chain within_comps begin
    @rtransform!(:area = (1 + :offset) * :occ)
    @rtransform!(:area = :area > 1.0 ? 1.0 : :area)
    @rtransform!(:area_per_site = :area / 500 * 100)
    @rtransform!(:site_per_area = (500 / 10_000) / :area)
end

# Visualize area
begin
    u = @rsubset(within_comps, :offset >= -0.5, :offset <= 0.5)
    @rtransform!(u, :side = first(split(:variable, "-")))
    f =
        data(u) *
        mapping(
            :value => "efficiency difference",
            :area;
            group=:sim => nonnumeric,
            color=:offset,
            row=:side,
        ) *
        visual(ScatterLines; colorscale=(x -> x + 0.5), markersize=3, linewidth=0.5)
    vline = mapping([0.0]) * visual(VLines; linestyle=:dash, color=:grey)
    fg = draw(f + vline, scales(; Color=(; colormap=:broc)))
    save(plotsdir("ranges_area_scatterlines.png"), fg)
    fg
end

# Area per site
begin
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(within_comps, :offset in set)
    rev = true

    # Random seed for jitter
    Random.seed!(42)

    # Main panel
    d1 = @rsubset(d, :set == "ranges")
    sortedcomps = unique(d1.variable)
    sortedoffsets = [
        o > 0.0 ? "+$(round(Int, 100o))" : "$(round(Int, 100o))" for o in unique(d1.offset)
    ]
    m = mapping(
        :variable =>
            renamer(sortedcomps .=> sortedoffsets) => "Range estimation difference (%)",
        :value => "Efficiency compared to True Range";
        color=:area_per_site,
    )
    rains = visual(
        RainClouds;
        markersize=6,
        jitter_width=0.30,
        plot_boxplots=false,
        clouds=nothing,
        orientation=:horizontal,
    )
    vline = mapping([0.0]) * visual(VLines; linestyle=:dash)
    hline =
        mapping([length(unique(d1.variable)) / 2 + 0.5]) *
        visual(HLines; linestyle=:solid, color=:lightgrey)
    scl = scales(; Color=(; colormap=:cividis))
    fg1 = draw(data(d1) * m * rains + vline + hline, scl; axis=(; xreversed=rev))
    # Figure
    save(plotsdir("ranges_area_scatter_area_per_site.png"), fg1)
    fg1
end

## Check the distribution of range efficiencies

# Sort values
sorted = unique(filter(:sim => ==(1), effs_estimations).variable)

# Median
med = median(filter(:variable => startswith("True"), effs_estimations).eff)

# Palette
palette = Dict()
for v in sorted
    if startswith("Over")(v)
        col = "1"
    elseif startswith("Under")(v)
        col = "2"
    else
        col = "3"
    end
    palette[v] = col
end

# Plot all efficiencies
begin
    Random.seed!(42) # for jitter
    eff = :eff => "efficiency"
    layer =
        mapping(
            :variable => sorter(sorted) => "", eff; color=:variable => (x -> palette[x])
        ) * visual(
            RainClouds;
            markersize=5,
            jitter_width=0.1,
            plot_boxplots=false,
            orientation=:horizontal,
        )
    vline = mapping([med]) * visual(VLines; linestyle=:dash)
    f1 = data(effs_estimations) * layer + vline
    f = Figure(; size=(700, 800))
    # sc = scales(; Color=(; palette=palette))
    ylog2 = (; axis=(;))
    fg1 = draw!(f[1, 1], f1; ylog2...)
    pad = (-45, 0, 20, 0)
    Label(f[1, 1, Top()], "Range estimations"; halign=:left, font=:bold, padding=pad)
    save(plotsdir("supp", "ranges_efficiencies_all.png"), f)
    f
end

# Plot subset
begin
    set = [
        "Under-0.50",
        "Under-0.40",
        "Under-0.30",
        "Under-0.20",
        "Under-0.10",
        "True-0.00",
        "Over-0.10",
        "Over-0.20",
        "Over-0.30",
        "Over-0.40",
        "Over-0.50",
    ]
    res = filter(:variable => in(set), effs_estimations)
    Random.seed!(42) # for jitter
    eff = :eff => "efficiency"
    layer =
        mapping(
            :variable => sorter(sorted) => "", eff; color=:variable => (x -> palette[x])
        ) * visual(
            RainClouds;
            markersize=5,
            jitter_width=0.1,
            # violin_limits=extrema,
            plot_boxplots=false,
            orientation=:horizontal,
            gap=0.01,
        )
    vline = mapping([med]) * visual(VLines; linestyle=:dash)
    f1 = data(res) * layer + vline
    f = Figure(;)
    # sc = scales(; Color=(; palette=palette))
    ylog2 = (; axis=(;))
    fg1 = draw!(f[1, 1], f1; ylog2...)
    pad = (-45, 0, 20, 0)
    Label(f[1, 1, Top()], "Range estimations"; halign=:left, font=:bold, padding=pad)
    save(plotsdir("supp", "ranges_efficiencies.png"), f)
    f
end

## Band figures for sign only (no overlap)

begin
    # Figure options
    f = Figure(; size=(900, 450))
    rev = false

    # Select results for bands
    res_bands = @rsubset(within_bands, :offset >= -0.5, :offset <= 0.5)
    var = :offset

    # Bands
    p1 = f[1:5, 1:3]
    ax = make_bands_ax!(p1; rev=false)
    Legend(f[:, end + 1], ax, "90% Percentile range"; framevisible=false)
    # ax = Axis(p1)

    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(within_comps, :offset in set)
    u = @rsubset(unique_comps, :offset in set, :countmeasure in ["sign_pos", "sign_neg"])

    # Comparison axis
    Random.seed!(42)
    col1 = Makie.wong_colors()[1]
    col2 = Makie.wong_colors()[2]
    colfunc(x) = [v < 0 ? col1 : col2 for v in x]
    scatter!(
        ax,
        [o + 0.02 * (rand() - 0.5) for o in d.offset],
        d.value;
        color=colfunc(d.value),
        markersize=5,
        alpha=0.7,
    )

    # Set axis limits
    ymin = minimum(res_bands.low)
    ymax = maximum(res_bands.upp)
    yoff = 0.1maximum(abs.([ymin, ymax]))
    ylims!(ax, minimum(res_bands.low) - yoff, maximum(res_bands.upp) + yoff)

    # Summary panel
    d3 = @rsubset(u, :set == "ranges")
    p2 = f[end + 1, 1:(end - 1)]
    ax0 = Axis(p2; yticks=([0.5, 1.5], ["Negative", "Positive"]))
    m34 = mapping(
        :offset,
        [1];
        stack=:countmeasure => sorter(["sign_neg", "sign_pos"]),
        color=:countmeasure => sorter(["sign_neg", "sign_pos"]),
        bar_labels=:label => verbatim,
    )
    v3 = visual(
        BarPlot;
        direction=:y,
        label_position=:center,
        label_color=:white,
        label_font=:bold,
        label_size=14,
        alpha=0.85,
    )
    draw!(ax0, data(d3) * m34 * v3)
    hidexdecorations!(ax0;)
    hideydecorations!(ax0; ticklabels=false, ticks=false)
    hidespines!(ax0)

    # Align axes
    linkxaxes!(ax, ax0)

    # Add labels
    pad = (-65, 0, 10, 0)
    Label(
        p1[1, 1, Top()],
        "Efficiency comparison between range estimations";
        halign=:left,
        font=:bold,
        padding=pad,
    )
    Label(p2[1, 1, Top()], "Sign summary"; font=:bold, halign=:left, padding=pad)

    # Figure
    save(plotsdir("ranges_efficiency_one.png"), f)
    f
end
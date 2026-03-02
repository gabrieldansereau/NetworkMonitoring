# using DrWatson
# @quickactivate :NetworkMonitoring

using Printf
include("include.jl") # see note regarding why we cannot use the module

# Load data
effs_estimations = CSV.read(datadir("efficiency_estimations.csv"), DataFrame)

# Rename
@transform!(
    effs_estimations,
    :variable = [
        @sprintf("%s-%.2f", s[1], parse(Float64, s[2])) for s in split.(:variable, "-")
    ]
)

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
    f = Figure(; size=(700, 700))
    # sc = scales(; Color=(; palette=palette))
    ylog2 = (; axis=(;))
    fg1 = draw!(f[1, 1], f1; ylog2...)
    pad = (-45, 0, 20, 0)
    Label(f[1, 1, Top()], "Range estimations"; halign=:left, font=:bold, padding=pad)
    f
end
save(plotsdir("ranges_efficiencies_all.png"), f)

# Plot subset
begin
    set = [
        "Under-0.30",
        "Under-0.20",
        "Under-0.10",
        "True-0.00",
        "Over-0.10",
        "Over-0.20",
        "Over-0.30",
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
            plot_boxplots=false,
            orientation=:horizontal,
        )
    vline = mapping([med]) * visual(VLines; linestyle=:dash)
    f1 = data(res) * layer + vline
    f = Figure(;)
    # sc = scales(; Color=(; palette=palette))
    ylog2 = (; axis=(;))
    fg1 = draw!(f[1, 1], f1; ylog2...)
    pad = (-45, 0, 20, 0)
    Label(f[1, 1, Top()], "Range estimations"; halign=:left, font=:bold, padding=pad)
    f
end
save(plotsdir("ranges_efficiencies.png"), f)

## Within-simulation comparison

# Separate results per simulation
within_combined_dif = comparewithin(
    effs_estimations,
    ["True-0.00", "Under-0.10", "Over-0.10"];
    labels=Dict("True-0.00" => "True"),
    f=(n, n2) -> efficiency_difference(n, n2; k=10_000),
)

# Count positive and negative comparisons per set and variable (across simulations/replicates)
unique_df = @chain within_combined_dif begin
    @groupby(:set, :variable)
    @combine(:count_pos = count(>(0), :value), :count_neg = count(<=(0), :value))
    stack([:count_pos, :count_neg]; variable_name=:countmeasure, value_name=:count)
end

# Visualize
sortedcomps = unique(unique_df.variable)
begin
    m = mapping(:variable, :value => "Δefficiency"; color=:value => (x -> x >= 0.0))
    rains = visual(
        RainClouds;
        markersize=7,
        jitter_width=0.15,
        plot_boxplots=false,
        clouds=nothing,
        orientation=:horizontal,
    )
    vline = mapping([0.0]) * visual(VLines; linestyle=:dash)
end
let d = within_combined_dif, u = unique_df
    Random.seed!(42)

    # Figure & grid
    f = Figure(; size=(600, 300))
    g1 = GridLayout(f[1:6, 1:3])
    g3 = GridLayout(f[1:6, end + 1])

    # Main panels
    d1 = @rsubset(d, :set == "ranges")
    m = mapping(
        :variable => sorter(sortedcomps) => "",
        :value => "Efficiency difference";
        color=:value => (x -> x <= 0.0),
    )
    # xlog2f = vs -> [rich("2", superscript("$(v)")) for v in vs]
    # xlog2 = (; axis=(; xtickformat=xlog2f))
    # xaxis = (; axis=(; xticks=0:2:10))
    fg1 = draw!(g1, data(d1) * m * rains + vline;)
    pad = (-100, 0, 10, 0)
    Label(
        g1[1, 1, Top()],
        "Range estimation comparisons";
        halign=:left,
        font=:bold,
        padding=pad,
    )

    # Summary panels
    d3 = @rsubset(u, :set == "ranges")
    ax3 = Axis(g3[1, 1])
    m34 = mapping(
        :variable => sorter(sortedcomps),
        [1];
        color=:countmeasure => sorter(["count_pos", "count_neg"]),
        stack=:countmeasure => sorter(["count_neg", "count_pos"]),
    )
    v3 = visual(
        BarPlot;
        direction=:x,
        bar_labels=["$v" for v in d3.count],
        label_position=:center,
        label_color=:white,
    )
    draw!(ax3, data(d3) * m34 * v3)

    hidedecorations!(ax3)
    hidespines!(ax3)
    pad = (0, 0, 10, 0)
    Label(g3[1, 1, Top()], "Frequency"; font=:bold, padding=pad)

    f
end
save(plotsdir("ranges_comparison.png"), current_figure())
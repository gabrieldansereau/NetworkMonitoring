# using DrWatson
# @quickactivate :NetworkMonitoring

using Printf
include("include.jl") # see note regarding why we cannot use the module

# Load data
effs_estimations = CSV.read(datadir("efficiency_estimations.csv"), DataFrame)

# Rename with same number of digits
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

## Within-simulation comparison

# Select the set
set = ["Under-0.20", "Under-0.10", "True-0.00", "Over-0.10", "Over-0.20"]

# Compare
within_combined_dif2 = comparewithin(
    select(effs_estimations, Not(["eff_low", "eff_upp"])),
    set;
    to="True-0.00",
    labels=Dict("True-0.00" => "True"),
    f=(n, n2) -> n - n2,
)
flipthatcomp!(within_combined_dif2, unique(within_combined_dif2.variable))
@rtransform!(
    within_combined_dif2, :variable = replace(:variable, "Δ" => "", "True" => "", "_" => "")
)

# Calculate intervals
effs_intervals = @chain effs_estimations begin
    @rsubset(:variable in set)
    @rtransform(:min = :eff_low - :eff, :max = :eff_upp - :eff)
    @select(:sim, :set, :variable, :min, :max)
end
effs_intervals_true = @rsubset(effs_intervals, :variable == "True-0.00")
true_min = Dict(r.sim => r.min for r in eachrow(effs_intervals_true))
true_max = Dict(r.sim => r.max for r in eachrow(effs_intervals_true))

# Check overlap
tmp = leftjoin(within_combined_dif2, effs_intervals; on=[:sim, :set, :variable])
newtmp = @chain tmp begin
    @select(Not(:occ))
    @rtransform(:true_min = true_min[:sim], :true_max = true_max[:sim])
    @rtransform(
        :overlap_neg = (:value + abs(:max) + abs(:true_min)) >= 0,
        :overlap_pos = (:value - abs(:min) - abs(:true_max)) <= 0,
    )
    @rtransform(:overlap = :value < 0 ? :overlap_neg : :overlap_pos)
end

# How'z it?
@rsubset(newtmp, :overlap == true)
@rsubset(newtmp, :overlap == true, abs(:value) > 1000)
@rsubset(newtmp, :overlap == false, abs(:value) < 100)
@chain begin
    newtmp
    @groupby(:variable)
    @combine(:pct = mean(:overlap))
end

# Count positive and negative comparisons per set and variable (across simulations/replicates)
unique_df2 = @chain newtmp begin
    @groupby(:set, :variable)
    @combine(
        :count_pos = count(>(0), :value) / length(:value),
        :count_neg = count(<=(0), :value) / length(:value)
    )
    stack([:count_pos, :count_neg]; variable_name=:countmeasure, value_name=:count)
end

# Visualize
begin
    d = newtmp
    u = unique_df2
    sortedcomps = unique(u.variable)
    rev = true

    # Figure & grid
    if nrow(u) <= 12
        f = Figure(; size=(750, 300))
    else
        f = Figure(; size=(600, 800))
    end

    function make_comps_ax!(f; d=d, u=u, sortedcomps=sortedcomps, l1, rev=rev)
        # Random seed for jitter
        Random.seed!(42)

        # Define grid layouts
        g1 = GridLayout(f[1:6, 1:3])
        g3 = GridLayout(f[1:6, end + 1])

        # Main panel
        d1 = @rsubset(d, :set == "ranges")
        m = mapping(
            :variable =>
                renamer([s => replace(s, "ΔTrue_" => "") for s in sortedcomps]) => "",
            :value => "Efficiency compared to True Range";
            # color=:value => (x -> x > 0.0),
            color=:overlap,
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
            mapping([(nrow(u) / 4) + 0.5]) *
            visual(HLines; linestyle=:solid, color=:lightgrey)
        col1 = Makie.wong_colors()[3]
        col2 = Makie.wong_colors()[4]
        fg1 = draw!(
            g1,
            data(d1) * m * rains + vline + hline,
            scales(; Color=(; palette=[col1, col2]));
            axis=(; xreversed=rev),
        )
        pad = (-80, 0, 10, 0)
        Label(g1[1, 1, Top()], l1; halign=:left, font=:bold, padding=pad)
        legend!(g1[1,2], fg1)

        # Add overlap
        pct = round(Int, 100 * mean(d1.overlap))
        Label(
            g1[1, 1, TopRight()],
            "Overlap: $pct %";
            tellheight=false,
            tellwidth=false,
            padding=pad,
        )
        for (i, v) in enumerate(sortedcomps)
            _df = @rsubset(d1, :variable == v)
            pct_v = round(Int, 100 * mean(_df.overlap))
            text!(
                -3500,
                i;
                text="$pct_v %",
                font=:bold,
                align=(:center, :baseline),
                color=pct_v <= 50 ? col1 : col2,
            )
        end

        # Summary panels
        d3 = @rsubset(u, :set == "ranges")
        ax3 = Axis(
            g3[1, 1];
            xreversed=!rev, # rev needs to be opposite somehow
            xticks=([0.5, 1.5], ["Negative", "Positive"]),
        )
        m34 = mapping(
            :variable => sorter(sortedcomps),
            [1];
            stack=:countmeasure => sorter(["count_neg", "count_pos"]),
            color=:countmeasure => sorter(["count_pos", "count_neg"]),
        )
        v3 = visual(
            BarPlot;
            direction=:x,
            bar_labels=["$(round(Int, v*100)) %" for v in d3.count],
            label_position=:center,
            label_color=:white,
            label_font=:bold,
            label_size=16,
            alpha=0.85,
        )
        draw!(ax3, data(d3) * m34 * v3)
        hideydecorations!(ax3)
        hidexdecorations!(ax3; ticklabels=false, ticks=false)
        hidespines!(ax3)

        # Align axes
        ax1 = g1.content[1].content
        linkyaxes!(ax1, ax3)

        # Add label
        pad = (0, 0, 10, 0)
        Label(g3[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)

        return (g1, g3)
    end
    make_comps_ax!(f; l1="Range estimation comparisons")
    save(plotsdir("ranges_confidence_interval.png"), current_figure())
    f
end

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
    f
end
save(plotsdir("ranges_efficiencies.png"), f)

## Within-simulation comparison

# Separate results per simulation
set = [
    # "Under-0.30",
    "Under-0.20",
    "Under-0.10",
    "True-0.00",
    "Over-0.10",
    "Over-0.20",
    # "Over-0.30",
]
select!(effs_estimations, Not(:eff_low, :eff_upp))
within_combined_dif2 = comparewithin(
    effs_estimations,
    set;
    to="True-0.00",
    labels=Dict("True-0.00" => "True"),
    # f=(n, n2) -> efficiency_difference(n, n2; k=10_000),
    f=(n, n2) -> n - n2,
)
flipthatcomp!(within_combined_dif2, unique(within_combined_dif2.variable))
@rtransform!(
    within_combined_dif2, :variable = replace(:variable, "Δ" => "", "True" => "", "_" => "")
)

# Count positive and negative comparisons per set and variable (across simulations/replicates)
unique_df2 = @chain within_combined_dif2 begin
    @groupby(:set, :variable)
    @combine(
        :count_pos = count(>=(0), :value) / length(:value),
        :count_neg = count(<(0), :value) / length(:value)
    )
    stack([:count_pos, :count_neg]; variable_name=:countmeasure, value_name=:count)
    @rtransform(:label = "$(round(Int, :count *100)) %")
    @rtransform(:label = (:count > 0.0 && :label == "0 %") ? "< 1 %" : :label)
end

# Visualize
begin
    d = within_combined_dif2
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
            color=:value => (x -> x > 0.0),
        )
        rains = visual(
            RainClouds;
            markersize=7,
            jitter_width=0.15,
            plot_boxplots=false,
            clouds=nothing,
            orientation=:horizontal,
        )
        vline = mapping([0.0]) * visual(VLines; linestyle=:dash)
        hline =
            mapping([(nrow(u) / 4) + 0.5]) *
            visual(HLines; linestyle=:solid, color=:lightgrey)
        fg1 = draw!(g1, data(d1) * m * rains + vline + hline; axis=(; xreversed=rev))
        pad = (-80, 0, 10, 0)
        Label(g1[1, 1, Top()], l1; halign=:left, font=:bold, padding=pad)

        # Summary panels
        d3 = @rsubset(u, :set == "ranges")
        ax3 = Axis(
            g3[1, 1];
            # xreversed=!rev, # rev needs to be opposite somehow
            xticks=([0.5, 1.5], ["Positive", "Negative"]),
        )
        m34 = mapping(
            :variable => sorter(sortedcomps),
            [1];
            stack=:countmeasure => sorter(["count_pos", "count_neg"]),
            color=:countmeasure => sorter(["count_neg", "count_pos"]),
            bar_labels=:label => verbatim,
        )
        v3 = visual(
            BarPlot;
            direction=:x,
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
        # Option 1
        # ax3.xticklabelpad = 7.0 # aligned with ax1
        # ax3.xticklabelpad = -3.0 # reasonable
        # Option 2
        # ax3.xaxisposition = :top
        # hidexdecorations!(ax3; ticklabels=false)
        # ax3.xticklabelpad = -5.0

        # Add label
        pad = (0, 0, 10, 0)
        Label(g3[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)

        return (g1, g3)
    end
    make_comps_ax!(f; l1="Range estimation comparisons")
    f
end
save(plotsdir("ranges_comparison.png"), current_figure())

## Plot bands (on the run!)

# Complete set of comparison
set_all = reverse(unique(effs_estimations.variable))
within_combined_all = comparewithin(
    effs_estimations,
    set_all;
    to="True-0.00",
    labels=Dict("True-0.00" => "True"),
    # f=(n, n2) -> efficiency_difference(n, n2; k=10_000),
    f=(n, n2) -> n - n2,
)
flipthatcomp!(within_combined_all, unique(within_combined_all.variable))
@rtransform!(
    within_combined_all, :variable = replace(:variable, "Δ" => "", "True" => "", "_" => "")
)

# Extract quantile range
within_bands = @chain within_combined_all begin
    groupby([:set, :variable])
    @combine(
        :low = quantile(:value, 0.05), :med = median(:value), :upp = quantile(:value, 0.95),
    )
    @select(:set, :variable, :offset = :variable, All())
    @rtransform(:offset = parse(Float64, replace(:offset, "Over-" => "", "Under" => "")))
end

# Add an Entry for offset of zero
push!(
    within_bands,
    (; set="ranges", variable="True-0.0", offset=0.0, low=0.0, med=0.0, upp=0.0),
)
sort!(within_bands, :offset)

# Restrict offset limits
@rsubset!(within_bands, :offset >= -0.24, :offset <= 0.24)

# Bands
tickdict = Dict(within_bands.offset .=> within_bands.variable)
tickdict[0.0] = "0.0"
fig_types = begin
    res = within_bands
    var = :offset
    rev = true

    fig = Figure(; size=(750, 350))
    function make_bands_ax!(f; res=res, var=var, rev=rev)
        t1, t2 = extrema(res[:, var])
        ax = Axis(
            f;
            xlabel="Offset",
            ylabel="Efficiency difference with True Range",
            xticks=ceil(t1; digits=1):0.1:floor(t2; digits=1),
            xtickformat=x -> [tickdict[tick] for tick in x],
            yreversed=rev,
            # xticklabelrotation=pi / 8,
        )
        if !rev
            limits!(ax, (nothing, nothing), (-1500, 1500))
        end
        x = res.offset
        low = res.low
        med = res.med
        upp = res.upp

        # Single band color option
        # col = Makie.wong_colors()[1]
        # band!(ax, x, low, upp; alpha=0.4, label=lab, color=col)
        # lines!(ax, x, low; linewidth=0.5, alpha=0.5, color=col)
        # lines!(ax, x, upp; linewidth=0.5, alpha=0.5, color=col)
        # lines!(ax, x, med; label="Median", color=col)

        # Two band color option
        col1 = Makie.wong_colors()[1]
        col2 = Makie.wong_colors()[2]
        lowpos = [l < 0 ? 0.0 : l for l in low]
        lowneg = [l < 0 ? l : 0.0 for l in low]
        upppos = [l < 0 ? 0.0 : l for l in upp]
        uppneg = [l < 0 ? l : 0.0 for l in upp]
        band!(ax, x, lowneg, uppneg; alpha=0.6, label="Lower efficiency", color=col1)
        band!(ax, x, lowpos, upppos; alpha=0.6, label="Higher efficiency", color=col2)
        cf(x) = [v < 0 ? col1 : col2 for v in x]
        lines!(ax, x, low; linewidth=0.5, alpha=0.5, color=cf(low))
        lines!(ax, x, upp; linewidth=0.5, alpha=0.5, color=cf(upp))
        lines!(ax, x, med; label="Median", color=col1)

        # Common options
        vlines!(ax, 0.0; linestyle=:solid, color=:lightgrey)
        hlines!(ax, 0.0; linestyle=:dash, color=:black)

        # Legend
        # axislegend(; position=:rt, merge=true)
        # axislegend(ax, "90% Interpercentile range"; position=:lt, merge=true)

        return ax
    end
    ax = make_bands_ax!(fig[1, 1])
    Legend(fig[1, 2], ax, "90% Percentile range")
    fig
end
save(plotsdir("ranges_bands.png"), current_figure())

# Combine subpanels
begin
    rev = true
    f = Figure(; size=(750, 650))
    g1, g3 = make_comps_ax!(
        f;
        d=within_combined_dif2,
        u=unique_df2,
        sortedcomps=unique(u.variable),
        rev=rev,
        l1="A) Efficiency comparison between range estimations",
    )
    ax2 = make_bands_ax!(
        f[(end + 1):(end + 10), 1:(end - 1)]; res=within_bands, var=:offset, rev=rev
    )
    Legend(f[7:end, end], ax2, "90% Percentile range"; framevisible=false)
    Label(
        f[7, 1, Top()],
        "B) Percentile range of efficiency differences";
        halign=:left,
        font=:bold,
        padding=(-80, 0, 10, 0),
    )
    save(plotsdir("ranges_combined.png"), f)
    f
end

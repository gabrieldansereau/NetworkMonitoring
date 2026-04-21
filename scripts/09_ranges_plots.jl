# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module

# Load data
effs_estimations = CSV.read(datadir("efficiency_estimations.csv"), DataFrame)

## Fill-in all possible offset values for simulations with missing results

# Check the proportion of simulations with results
@chain effs_estimations begin
    groupby([:variable])
    combine(nrow)
    @transform(:prop = :nrow ./ (maximum(:nrow)))
    @rsubset(:prop < 1.0)
end

# Get all combinations
effs_estimations_all = allcombinations(
    DataFrame;
    sim=unique(effs_estimations.sim),
    set=unique(effs_estimations.set),
    variable=unique(effs_estimations.variable),
)
leftjoin!(effs_estimations_all, effs_estimations; on=[:sim, :set, :variable])

# Check the missing results
@chain effs_estimations_all begin
    @rsubset(ismissing(:eff))
end

# Replace missing values by maximum
sims_missing = unique(@rsubset(effs_estimations_all, ismissing(:eff)).sim)
for sim in sims_missing
    # Select simulation
    effs_sim = @view effs_estimations_all[effs_estimations_all.sim .== sim, :]
    # Extract efficiency from highest non-missing offset
    _, max_i = findmax(skipmissing(effs_sim.offset))
    max_row = effs_sim[max_i, :]
    # Replace values
    for r in eachrow(effs_sim)
        if ismissing(r.eff)
            r.offset = parse(Float64, replace(r.variable, "Over-" => "", "Under" => ""))
            r.eff = max_row.eff[1]
            r.eff_low = max_row.eff_low[1]
            r.eff_upp = max_row.eff_upp[1]
            r.occ = max_row.occ[1]
            r.deg = max_row.deg
            r.pmax = max_row.pmax
        end
    end
end
disallowmissing!(effs_estimations_all)

# Use dataset with all possible values
use_all = false
if use_all
    effs_estimations = effs_estimations_all
end

## Within-simulation comparison

# Complete set of comparisons
set_all = reverse(unique(effs_estimations.variable))
within_comps = comparewithin(
    effs_estimations,
    set_all;
    to="True-0.00",
    labels=Dict("True-0.00" => "True"),
    f=(n, n2) -> n - n2,
)
flipthatcomp!(within_comps, unique(within_comps.variable))
@rtransform!(
    within_comps, :variable = replace(:variable, "Δ" => "", "True" => "", "_" => "")
)

# Bring back offset & update variable label
@rtransform!(
    within_comps, :offset = parse(Float64, replace(:variable, "Over-" => "", "Under" => ""))
)
select!(within_comps, :sim, :variable, :offset, :value, :overlap, Not(:set), All(), :set)

# Count positive and negative comparisons per set and variable (across simulations/replicates)
unique_comps = @chain within_comps begin
    @groupby :set :variable :offset
    @combine(
        :n = length(:value),
        :sign_pos = count(>=(0), :value) / length(:value),
        :sign_neg = count(<(0), :value) / length(:value),
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

# Extract quantile range for all offset
within_bands = @chain within_comps begin
    groupby([:set, :variable, :offset])
    @combine(
        :low = quantile(:value, 0.05), :med = median(:value), :upp = quantile(:value, 0.95),
    )
end

# Add an Entry for offset of zero
push!(
    within_bands,
    (; set="ranges", variable="True-0.0", offset=0.0, low=0.0, med=0.0, upp=0.0),
)
sort!(within_bands, :offset)

## Plot comparisons and bands

# Combined figure
begin
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(within_comps, :offset in set)
    u = @rsubset(unique_comps, :offset in set, :countmeasure in ["sign_pos", "sign_neg"])

    # Select results for bands
    res_bands = @rsubset(within_bands, :offset >= -0.5, :offset <= 0.5)
    var = :offset

    # Figure options
    f = Figure(; size=(800, 750))
    rev = false

    # Comparison panel
    function make_comps_ax!(g1, g3; d=d, u=u, l1, rev=rev)
        # Random seed for jitter
        Random.seed!(42)

        # Sorted sets for comparison
        sort!(u, :offset)
        sortedcomps = unique(u.variable)
        sortedoffsets = [
            o > 0.0 ? "+$(round(Int, 100o))" : "$(round(Int, 100o))" for
            o in unique(u.offset)
        ]
        # Main panel
        d1 = @rsubset(d, :set == "ranges")
        m = mapping(
            :variable =>
                renamer(sortedcomps .=> sortedoffsets) => "Range estimation difference (%)",
            :value => "Efficiency compared to True Range";
            color=:value => (x -> x > 0.0),
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
            mapping([length(unique(u.variable)) / 2 + 0.5]) *
            visual(HLines; linestyle=:solid, color=:lightgrey)
        fg1 = draw!(g1, data(d1) * m * rains + vline + hline; axis=(; xreversed=rev))
        pad = (-65, 0, 10, 0)
        Label(g1[1, 1, Top()], l1; halign=:left, font=:bold, padding=pad)
        Label(
            g1[1, 1, TopLeft()],
            "Over ⬆️";
            padding=(0, 2, -40, 0),
            halign=:right,
            fontsize=12,
        )
        Label(
            g1[1, 1, BottomLeft()],
            "Under ⬇️";
            padding=(0, 2, 0, -50),
            halign=:right,
            fontsize=12,
        )

        # Summary panels
        d3 = @rsubset(u, :set == "ranges")
        ax3 = Axis(g3[1, 1]; xticks=([0.5, 1.5], ["Negative", "Positive"]))
        m34 = mapping(
            :variable => sorter(sortedcomps),
            [1];
            stack=:countmeasure => sorter(["sign_neg", "sign_pos"]),
            color=:countmeasure => sorter(["sign_neg", "sign_pos"]),
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

        # Add labels
        pad = (0, 0, 10, 0)
        Label(g3[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)

        return (g1, g3)
    end

    # Bands panel
    function make_bands_ax!(
        f; res=res_bands, var=var, rev=rev, colour=true, arrowlabels=true
    )
        # Values
        x = res[:, var]
        low = res.low
        med = res.med
        upp = res.upp

        # Axis
        t1, t2 = extrema(x)
        tickdict = Dict(within_bands.offset .=> string.(within_bands.offset))
        tickdict[0.0] = "0.0"
        ax = Axis(
            f;
            xlabel="Range estimation difference (%)",
            ylabel="Efficiency compared to True Range",
            xticks=ceil(t1; digits=1):0.1:floor(t2; digits=1),
            xtickformat=values ->
                [v > 0.0 ? "+$(Int(100*v))" : "$(Int(100*v))" for v in values],
            yreversed=rev,
        )
        if !rev
            # limits!(ax, (nothing, nothing), (-1500, 1500))
        end

        # Two band color option
        if colour
            col1 = Makie.wong_colors()[1]
            col2 = Makie.wong_colors()[2]
            lab1 = "Lower efficiency"
            lab2 = "Higher efficiency"
        else
            col1 = :grey
            col2 = :grey
            lab1 = "90% Percentile range"
            lab2 = lab1
        end
        lowpos = [l < 0 ? 0.0 : l for l in low]
        lowneg = [l < 0 ? l : 0.0 for l in low]
        upppos = [l < 0 ? 0.0 : l for l in upp]
        uppneg = [l < 0 ? l : 0.0 for l in upp]
        band!(ax, x, lowneg, uppneg; alpha=0.6, label=lab1, color=col1)
        band!(ax, x, lowpos, upppos; alpha=0.6, label=lab2, color=col2)
        cf(x) = [v < 0 ? col1 : col2 for v in x]
        lines!(ax, x, low; linewidth=0.5, alpha=0.5, color=cf(low))
        lines!(ax, x, upp; linewidth=0.5, alpha=0.5, color=cf(upp))
        lines!(ax, x, med; label="Median", color=col1)

        # Common options
        vlines!(ax, 0.0; linestyle=:solid, color=:lightgrey)
        hlines!(ax, 0.0; linestyle=:dash, color=:black)

        # Add labels
        if arrowlabels
            Label(
                f[1, 1, BottomLeft()],
                "⬅️ Under";
                padding=(0, 0, -72, 0),
                halign=:left,
                fontsize=12,
                tellheight=false,
                tellwidth=false,
            )
            Label(
                f[1, 1, BottomRight()],
                "Over ➡️";
                padding=(0, 0, -72, 0),
                halign=:right,
                fontsize=12,
                tellheight=false,
                tellwidth=false,
            )
        end

        return ax
    end

    # Create figure
    g1 = GridLayout(f[1:6, 1:3])
    g3 = GridLayout(f[1:6, end + 1])
    g1, g3 = make_comps_ax!(
        g1, g3; rev=rev, l1="A) Efficiency comparison between range estimations"
    )
    ax2 = make_bands_ax!(f[(end + 1):(end + 5), 1:(end - 1)]; rev=rev)
    Legend(f[7:end, end], ax2, "90% Percentile range"; framevisible=false)
    Label(
        f[7, 1, Top()],
        "B) Percentile range of efficiency differences";
        halign=:left,
        font=:bold,
        padding=(-65, 0, 10, 0),
    )
    # Align Axis labels
    ax1 = content(g1[1, 1])
    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2]) + 2
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace
    # Save
    save(plotsdir("ranges_efficiency.png"), f)
    f
end

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

## Confidence intervals

# Extract quantile range for all offset
overlap_bands = rename(unique_comps, :count => :med)
select!(overlap_bands, Not(:label))

# Add an Entry for offset of zero
common = (;
    set="ranges", variable="True-0.0", offset=0.0, n=maximum(overlap_bands.n), prop_sim=1.0
)
push!(overlap_bands, (; common..., countmeasure="overlap", med=1.0))
push!(overlap_bands, (; common..., countmeasure="positive", med=0.0))
push!(overlap_bands, (; common..., countmeasure="negative", med=0.0))
sort!(overlap_bands, :offset)

# Add confidence interval for the proportion
@chain overlap_bands begin
    @rtransform!(:x = round(Int, :med * :n))
    @rtransform!(
        :low = confint(BinomialTest(:x, :n); level=0.90, method=(^(:wilson)))[1],
        :upp = confint(BinomialTest(:x, :n); level=0.90, method=(^(:wilson)))[2],
    )
    @select!(Not(:x))
end

# Visualize
begin
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(within_comps, :offset in set)
    u = @rsubset(
        unique_comps, :offset in set, :countmeasure in ["positive", "negative", "overlap"]
    )

    # Select results for bands
    res = @rsubset(overlap_bands, :offset in set)
    var = :offset

    # Figure options
    f = Figure(; size=(800, 750))
    rev = false

    # Overlap panel
    function make_overlap_ax!(g1, g3; d=d, u=u, l1, rev=rev, arrowlabels=true)
        # Random seed for jitter
        Random.seed!(42)

        # Sorted sets for comparison
        sort!(u, :offset)
        sortedcomps = unique(u.variable)
        sortedoffsets = [
            o > 0.0 ? "+$(round(Int, 100o))" : "$(round(Int, 100o))" for
            o in unique(u.offset)
        ]

        # Main panel
        d1 = @rsubset(d, :set == "ranges")
        m = mapping(
            :variable =>
                renamer(sortedcomps .=> sortedoffsets) => "Range estimation difference (%)",
            :value => "Efficiency compared to True Range";
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
            mapping([length(unique(u.variable)) / 2 + 0.5]) *
            visual(HLines; linestyle=:solid, color=:lightgrey)
        pal = [
            "negative" => Makie.wong_colors()[3],
            "overlap" => Makie.wong_colors()[4],
            "positive" => Makie.wong_colors()[2],
        ]
        scl = scales(; Color=(; palette=pal))
        fg1 = draw!(g1, data(d1) * m * rains + vline + hline, scl; axis=(; xreversed=rev))
        pad = (-65, 0, 10, 0)
        Label(g1[1, 1, Top()], l1; halign=:left, font=:bold, padding=pad)

        # Summary panels
        d3 = @rsubset(u, :set == "ranges")
        ax3 = Axis(g3[1, 1]; xticks=([0.5, 1.5, 2.5], ["Negative", "Overlap", "Positive"]))
        sortedmeasures = first.(pal)
        m34 = mapping(
            :variable => sorter(sortedcomps),
            [1];
            stack=:countmeasure => sorter(sortedmeasures),
            color=:countmeasure,
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
        draw!(ax3, data(d3) * m34 * v3, scl)
        hideydecorations!(ax3)
        hidexdecorations!(ax3; ticklabels=false, ticks=false)
        hidespines!(ax3)

        # Align axes
        ax1 = g1.content[1].content
        linkyaxes!(ax1, ax3)

        # Add labels
        pad = (0, 0, 10, 0)
        Label(g3[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)
        if arrowlabels
            Label(
                g1[1, 1, TopLeft()],
                "Over ⬆️";
                padding=(0, 2, -40, 0),
                halign=:right,
                fontsize=12,
            )
            Label(
                g1[1, 1, BottomLeft()],
                "Under ⬇️";
                padding=(0, 2, 0, -50),
                halign=:right,
                fontsize=12,
            )
        end

        return (g1, g3)
    end

    # Bands panel
    function make_overlap_bands!(f; res=res, var=var, rev=rev, title="", addlines=false)
        # Axis
        t1, t2 = extrema(set)
        ax = Axis(
            f[1, 1];
            xlabel="Range estimation difference (%)",
            ylabel="Proportion of simulations",
            xticks=ceil(t1; digits=1):0.1:floor(t2; digits=1),
            xtickformat=values ->
                [v > 0.0 ? "+$(Int(100*v))" : "$(Int(100*v))" for v in values],
            xreversed=rev,
            yticks=0.0:0.25:1.0,
            title=title,
        )
        # Band colors
        pal = Dict(
            "negative" => Makie.wong_colors()[3],
            "overlap" => Makie.wong_colors()[4],
            "positive" => Makie.wong_colors()[2],
        )
        # Bands
        for mes in ["overlap", "negative", "positive"]
            r = @rsubset(res, :countmeasure == mes)
            x = r[:, var]
            med = r.med
            low = r.low
            upp = r.upp
            lab = uppercasefirst(mes)
            col = pal[mes]
            band!(ax, x, low, upp; alpha=0.6, label=lab, color=col)
            lines!(ax, x, med; label=lab, color=col)
        end
        # Common options
        if addlines
            vlines!(ax, 0.0; linestyle=:solid, color=:grey)
            hlines!(ax, 0.0; linestyle=:dash, color=:black)
        end

        return ax
    end

    # Create figure
    g1 = GridLayout(f[1:6, 1:2])
    g3 = GridLayout(f[:, end + 1])
    g4l = 7:9
    make_overlap_ax!(g1, g3; l1="A) Efficiency comparison between range estimations")
    # ax2 = make_overlap_bands!(f[g4l, 1:2]; addlines=true)
    ax2 = make_overlap_bands!(
        f[g4l, 1]; res=@rsubset(res, :offset <= 0.0), rev=true, title="Underestimation"
    )
    ax3 = make_overlap_bands!(
        f[g4l, 2]; res=@rsubset(res, :offset >= 0.0), rev=false, title="Overestimation"
    )
    Legend(
        f[g4l, end], ax2, "Comparison sign"; framevisible=false, merge=true, tellwidth=false
    )
    Label(
        f[g4l, 1, Top()],
        "B) Change in proportion of comparisons";
        halign=:left,
        font=:bold,
        padding=(-65, 0, 30, 0),
    )
    # Align Axis labels
    ax1 = content(g1[1, 1])
    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2]) + 2
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace
    ax3.yticklabelspace = yspace
    # Save
    save(plotsdir("ranges_overlap.png"), current_figure())
    f
end

begin
    function make_comps_ax!(ax; d=d, res=res_bands)
        Random.seed!(42)
        scl = scales(; Color=(; palette=[k => v for (k, v) in pal]))
        sc = scatter!(
            ax,
            [o + 0.02 * (rand() - 0.5) for o in d.offset],
            d.value;
            color=[pal[v] for v in d.overlap],
            label=[uppercasefirst(v) => (; color=pal[v], markersize=8) for v in d.overlap],
            markersize=5,
            alpha=0.9,
        )
        # Set axis limits
        ymin = minimum(res.low)
        ymax = maximum(res.upp)
        yabsmax = maximum(abs.([ymin, ymax]))
        yoff = 0.1yabsmax
        ylims!(ax, ymin - yoff, ymax + yoff)
        return sc
    end
    let
        f = Figure()
        ax = Axis(f[1, 1])
        sc = make_comps_ax!(ax; d=d, res=res_bands)
        f
    end
end

begin
    function make_summary_ax!(gp; u=u)
        ax0 = Axis(gp; yticks=([0.5, 1.5, 2.5], ["Negative", "Overlap", "Positive"]))
        m34 = mapping(
            :offset,
            [1];
            stack=:countmeasure => sorter(["negative", "overlap", "positive"]),
            color=:countmeasure,
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
        draw!(ax0, data(u) * m34 * v3, scl)
        hidexdecorations!(ax0;)
        hideydecorations!(ax0; ticklabels=false, ticks=false)
        hidespines!(ax0)
        return ax0
    end
    let
        f = Figure(; size=(700, 200))
        g = GridLayout(f[:, :])
        ax = make_summary_ax!(g[:, :]; u=u)
        f
    end
end

begin
    # Figure options
    f = Figure(; size=(900, 450))
    rev = false

    # Select results for bands
    res_bands = @rsubset(within_bands, :offset >= -0.5, :offset <= 0.5)
    var = :offset
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(within_comps, :offset in set)
    u = @rsubset(
        unique_comps, :offset in set, :countmeasure in ["positive", "negative", "overlap"]
    )
    # Define colour palette
    pal = Dict(
        "negative" => Makie.wong_colors()[3],
        "overlap" => Makie.wong_colors()[4],
        "positive" => Makie.wong_colors()[2],
    )

    # GridLayout
    g1 = GridLayout(f[1:5, 1:3])
    g2 = GridLayout(f[1:5, 4])
    g3 = GridLayout(f[6:7, 1:3])

    # Bands
    ax1 = make_bands_ax!(g1[:, :]; res=res_bands, var=var, rev=false, colour=false)
    # Comparison axis
    sc1 = make_comps_ax!(ax1; d=d, res=res_bands)
    # Legend
    Legend(g2[:, :], ax1; framevisible=false, merge=true)
    # Summary panel
    ax3 = make_summary_ax!(g3[:, :]; u=u)
    # Align axes
    linkxaxes!(ax1, ax3)

    # Add labels
    pad = (-65, 0, 10, 0)
    Label(
        g1[1, 1, Top()],
        "Efficiency comparison between range estimations";
        halign=:left,
        font=:bold,
        padding=pad,
    )
    Label(g3[1, 1, Top()], "Sign summary"; font=:bold, halign=:left, padding=pad)

    # Figure
    save(plotsdir("ranges_overlap_one.png"), f)
    f
end

begin
    # Figure options
    f = Figure(; size=(900, 700))

    # GridLayout
    g1 = GridLayout(f[1:5, 1:3])
    g1l = GridLayout(f[1:5, 4])
    g2 = GridLayout(f[-1:0, 1:3])
    g3 = GridLayout(f[(end + 1):(end + 4), 1:3])
    g3l = GridLayout(f[(end - 3):end, 4])

    # Bands
    ax1 = make_bands_ax!(g1[:, :]; res=res_bands, var=var, rev=false, colour=false)
    # Comparison axis
    sc1 = make_comps_ax!(ax1; d=d, res=res_bands)
    # Legend
    Legend(g1l[:, :], ax1; framevisible=false, merge=true, halign=:left)
    # Summary panel
    ax2 = make_summary_ax!(g2[:, :]; u=u)
    vlines!(ax2, [0.0]; linestyle=:solid, color=:lightgrey)
    # ax2.xticksvisible=true
    # ax2.xticks=[-0.5:0.1:-0.1..., 0.1:0.1:0.5...]
    ax1.xticksmirrored = true
    # Overlap bands
    ax3 = make_overlap_bands!(g3[:, :]; res=overlap_bands, rev=false, title="")
    vlines!(ax3, [0.0]; linestyle=:solid, color=:lightgrey)
    Legend(g3l[:, :], ax3; framevisible=false, merge=true, tellwidth=false, halign=:left)

    # Align axes
    linkxaxes!(ax1, ax2, ax3)

    # Add labels
    pad = (-65, 0, 10, 0)
    Label(
        g2[1, 1, Top()],
        "a) Comparison of efficiencies between range estimations";
        halign=:left,
        font=:bold,
        padding=pad,
    )
    # Label(g2[1, 1, Top()], "Sign summary"; font=:bold, halign=:left, padding=pad)
    Label(
        g3[1, 1, Top()],
        "b) Proportion of simulations per comparison sign";
        font=:bold,
        halign=:left,
        padding=pad,
    )

    # Figure
    save(plotsdir("ranges_overlap_one_and_bands.png"), f)
    f
end

begin
    # Figure options
    f = Figure(; size=(900, 700))

    # GridLayout
    g1 = GridLayout(f[1:5, 1:3])
    g1l = GridLayout(f[1:5, 4])
    # g2 = GridLayout(f[-1:0, 1:3])
    g3 = GridLayout(f[(end + 1):(end + 4), 1:3])
    g3l = GridLayout(f[(end - 3):end, 4])

    # Bands
    ax1 = make_bands_ax!(g1[:, :]; res=res_bands, var=var, rev=false, colour=false)
    # Comparison axis
    sc1 = make_comps_ax!(ax1; d=d, res=res_bands)
    # Summary panel
    ax2 = make_summary_ax!(g2[:, :]; u=u)
    vlines!(ax2, [0.0]; linestyle=:solid, color=:lightgrey)
    # Overlap bands
    ax3 = make_overlap_bands!(g3[:, :]; res=overlap_bands, rev=false, title="")
    vlines!(ax3, [0.0]; linestyle=:solid, color=:lightgrey)

    # Legends
    leg_opt = (; framevisible=false, merge=true, halign=:left)
    Legend(g1l[:, :], ax1, "Comparison values"; leg_opt...)
    Legend(g3l[:, :], ax3, "Comparison sign"; leg_opt..., tellwidth=false)

    # Align axes
    linkxaxes!(ax1, ax2, ax3)

    # Add labels
    pad = (-65, 0, 10, 0)
    lab_opt = (; halign=:left, font=:bold, padding=pad)
    Label(
        g1[1, 1, Top()],
        "a) Comparison of efficiencies between range estimations";
        lab_opt...,
    )
    Label(g3[1, 1, Top()], "b) Proportion of simulations per comparison sign"; lab_opt...)

    # Figure
    save(plotsdir("ranges_overlap_minimal.png"), f)
    f
end

begin
    # Figure options
    f = Figure(; size=(900, 600))

    # GridLayout
    g1 = GridLayout(f[1:5, 1:4])
    g1a = GridLayout(g1[:, 1:2])
    g1b = GridLayout(g1[:, 3:4])
    g1l = GridLayout(g1[:, 5])
    g3 = GridLayout(f[(end + 1):(end + 4), :])
    g3a = GridLayout(g3[:, 1:2])
    g3b = GridLayout(g3[:, 3:4])
    g3l = GridLayout(g3[:, 5])

    # Bands
    res1a = @rsubset(res_bands, :offset <= 0.0)
    res1b = @rsubset(res_bands, :offset >= 0.0)
    ax1a = make_bands_ax!(
        g1a[:, :]; res=res1a, var=var, rev=false, colour=false, arrowlabels=false
    )
    ax1b = make_bands_ax!(
        g1b[:, :]; res=res1b, var=var, rev=false, colour=false, arrowlabels=false
    )
    hideydecorations!(ax1b; grid=false)
    # Comparison axis
    d1a = @rsubset(d, :offset <= 0.0)
    d1b = @rsubset(d, :offset >= 0.0)
    sc1a = make_comps_ax!(ax1a; d=d1a, res=res1a)
    sc1b = make_comps_ax!(ax1b; d=d1b, res=res1b)
    hidexdecorations!(ax1a; grid=false)
    hidexdecorations!(ax1b; grid=false)
    # Overlap bands
    res3a = @rsubset(overlap_bands, :offset <= 0.0)
    res3b = @rsubset(overlap_bands, :offset >= 0.0)
    ax3a = make_overlap_bands!(g3a[:, :]; res=res3a, rev=false)
    ax3b = make_overlap_bands!(g3b[:, :]; res=res3b, rev=false)
    vlines!(ax3a, [0.0]; linestyle=:solid, color=:lightgrey)
    vlines!(ax3b, [0.0]; linestyle=:solid, color=:lightgrey)
    hideydecorations!(ax3b; grid=false)

    # Legends
    leg_opt = (; framevisible=false, merge=true, halign=:left)
    Legend(g1l[:, :], ax1a, "Comparison values"; leg_opt..., tellwidth=false)
    Legend(g3l[:, :], ax3a, "Comparison sign"; leg_opt..., tellwidth=false)

    # Adjust axis limits
    ymin = minimum(res_bands.low)
    ymax = maximum(res_bands.upp)
    yabsmax = maximum(abs.([ymin, ymax]))
    yoff = 0.1yabsmax
    ylims!(ax1a, ymin - yoff, ymax + yoff)
    ylims!(ax1b, ymin - yoff, ymax + yoff)

    # Titles
    ax1a.title = "Underestimation"
    ax1b.title = "Overestimation"

    # Align axes
    # ax1a.xreversed=true
    # ax3a.xreversed=true
    linkxaxes!(ax1a, ax3a)
    linkxaxes!(ax1b, ax3b)
    linkyaxes!(ax1a, ax1b)
    linkyaxes!(ax3a, ax3b)

    # Add labels
    pad = (-65, 0, 30, 0)
    lab_opt = (; halign=:left, font=:bold, padding=pad)
    # Label(
    #     g1[1, 1, Top()],
    #     "a) Comparison of efficiencies between range estimations";
    #     lab_opt...,
    # )
    # Label(g3[1, 1, Top()], "b) Proportion of simulations per comparison sign"; lab_opt...)

    # Figure
    save(plotsdir("ranges_overlap_minimal_split.png"), f)
    f
end

## Mix overlap & efficiency sign

begin
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(within_comps, :offset in set)
    u = @rsubset(
        unique_comps, :offset in set, :countmeasure in ["positive", "negative", "overlap"]
    )

    # Select results for bands
    res_bands = @rsubset(within_bands, :offset >= -0.5, :offset <= 0.5)
    var = :offset

    # Figure options
    f = Figure(; size=(800, 750))
    rev = false

    # Create figure
    g1 = GridLayout(f[1:6, 1:2])
    g3 = GridLayout(f[:, end + 1])
    g1, g3 = make_overlap_ax!(
        g1, g3; d=d, u=u, rev=rev, l1="A) Efficiency comparison between range estimations"
    )
    ax2 = make_bands_ax!(
        f[(end + 1):(end + 5), 1:(end - 1)]; res=res_bands, var=var, rev=rev
    )
    Legend(f[7:end, end], ax2, "90% Percentile range"; framevisible=false, width=200)
    Label(
        f[7, 1, Top()],
        "B) Percentile range of efficiency differences";
        halign=:left,
        font=:bold,
        padding=(-65, 0, 10, 0),
    )
    # Align Axis labels
    ax1 = content(g1[1, 1])
    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2]) + 10
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace
    # Save
    save(plotsdir("ranges_mixed.png"), current_figure())
    f
end

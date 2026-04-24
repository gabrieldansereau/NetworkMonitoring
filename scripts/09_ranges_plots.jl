# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module

# Load data
pmax_opt = "n_at_pmax4"
effs_estimations = CSV.read(datadir("efficiency_estimations-$pmax_opt.csv"), DataFrame)

# Inverse lower and upper bounds
@chain effs_estimations begin
    @rename!(:eff_low1 = :eff_low, :eff_upp1 = :eff_upp)
    @rename!(:eff_low = :eff_upp1, :eff_upp = :eff_low1)
end

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
within_combined_all = comparewithin(
    select(effs_estimations, :sim, :set, :variable, :eff, :occ),
    set_all;
    to="True-0.00",
    labels=Dict("True-0.00" => "True"),
    f=(n, n2) -> n - n2,
)
flipthatcomp!(within_combined_all, unique(within_combined_all.variable))
@rtransform!(
    within_combined_all, :variable = replace(:variable, "Δ" => "", "True" => "", "_" => "")
)

# Bring back offset & update variable label
@rtransform!(
    within_combined_all,
    :offset = parse(Float64, replace(:variable, "Over-" => "", "Under" => ""))
)
select!(within_combined_all, Not(:value), :value)

# Count positive and negative comparisons per set and variable (across simulations/replicates)
unique_comps = @chain within_combined_all begin
    @groupby(:set, :variable, :offset)
    @combine(
        :prop_sim = length(:value),
        :count_pos = count(>=(0), :value) / length(:value),
        :count_neg = count(<(0), :value) / length(:value)
    )
    @transform(:prop_sim = :prop_sim ./ maximum(:prop_sim))
    stack([:count_pos, :count_neg]; variable_name=:countmeasure, value_name=:count)
    @rtransform(:label = "$(round(Int, :count *100)) %")
    @rtransform(:label = (:count > 0.0 && :label == "0 %") ? "< 1 %" : :label)
end

# Extract quantile range for all offset
within_bands = @chain within_combined_all begin
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

## Add sampling effort over area

#=

# Add area and sampling effort
@chain within_combined_all begin
    @rtransform!(:area = (1 + :offset) * :occ)
    @rtransform!(:area = :area > 1.0 ? 1.0 : :area)
    @rtransform!(:area_per_site = :area / 500 * 100)
    @rtransform!(:site_per_area = (500 / 10_000) / :area)
end

# Visualize area
begin
    u = @rsubset(within_combined_all, :offset >= -0.5, :offset <= 0.5)
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
    d = @rsubset(within_combined_all, :offset in set)
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

=#

## Plot comparisons and bands

# Combined figure
begin
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(within_combined_all, :offset in set)
    u = @rsubset(unique_comps, :offset in set)

    # Select results for bands
    res_bands = @rsubset(within_bands, :offset >= -0.5, :offset <= 0.5)
    var = :offset

    # Figure options
    f = Figure(; size=(800, 750))
    rev = true

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

        # Summary panels
        d3 = @rsubset(u, :set == "ranges")
        ax3 = Axis(g3[1, 1]; xticks=([0.5, 1.5], ["Positive", "Negative"]))
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

        # Add labels
        pad = (0, 0, 10, 0)
        Label(g3[1, 1, Top()], "Comparison sign"; font=:bold, padding=pad)
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

        return (g1, g3)
    end

    # Bands panel
    function make_bands_ax!(f; res=res_bands, var=var, rev=rev)
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
            ylabel="Efficiency difference with True Range",
            xticks=ceil(t1; digits=1):0.1:floor(t2; digits=1),
            xtickformat=values ->
                [v > 0.0 ? "+$(Int(100*v))" : "$(Int(100*v))" for v in values],
            yreversed=rev,
        )
        if !rev
            limits!(ax, (nothing, nothing), (-1500, 1500))
        end

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

        # Add labels
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
    Label(f[end:end, end, Bottom()], string(pmax_opt); font=:bold)
    # Align Axis labels
    ax1 = content(g1[1, 1])
    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2]) + 2
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace
    # Save
    save(plotsdir("ranges_efficiency.png"), f)
    f
end

## Confidence intervals

# Calculate intervals
effs_intervals_true = @rsubset(effs_estimations, :variable == "True-0.00")
true_low = Dict(r.sim => r.eff_low for r in eachrow(effs_intervals_true))
true_upp = Dict(r.sim => r.eff_upp for r in eachrow(effs_intervals_true))

# Check overlap
effs_overlap = @chain effs_estimations begin
    # Remove the True Range simulations
    @rsubset(:offset != 0.0)
    # Join back the comparison values with the actual efficiency values
    leftjoin(within_combined_all; on=[:sim, :set, :variable, :offset, :occ])
    @select(:sim, :variable, :offset, :value, :eff_low, :eff_upp, All())
    # Check the overlap with the interval for the True Range
    @rtransform(:true_low = true_low[:sim], :true_upp = true_upp[:sim])
    @rtransform(:overlap = :eff_low <= :true_upp && :true_low <= :eff_upp)
    # Assess the overlap sign
    @rtransform(
        :overlap_sign =
            :overlap == false ? (:value < 0 ? "negative" : "positive") : "overlap"
    )
    # Filter out columns used for intervals and order nicely
    @select(
        :sim, :variable, :offset, :value, :overlap, :overlap_sign, :deg, :pmax, :occ, :set
    )
end

# Count positive and negative comparisons per set and variable (across simulations/replicates)
unique_overlap = @chain effs_overlap begin
    @groupby(:set, :variable, :offset)
    @combine(
        :n = length(:value),
        :positive = count(==("positive"), :overlap_sign) / length(:overlap_sign),
        :negative = count(==("negative"), :overlap_sign) / length(:overlap_sign),
        :overlap = count(==("overlap"), :overlap_sign) / length(:overlap_sign),
    )
    stack([:positive, :negative, :overlap]; variable_name=:countmeasure, value_name=:count)
    @rtransform(:label = "$(round(Int, :count *100)) %")
    @rtransform(:label = (:label == "0 %" && :count > 0.0) ? "< 1 %" : :label)
    @rtransform(:label = (:label == "1 %" && :count < 1.0) ? "< 1 %" : :label)
end

# Extract quantile range for all offset
overlap_bands = rename(unique_overlap, :count => :med)
select!(overlap_bands, Not(:label))

# Add an Entry for offset of zero
common = (; set="ranges", variable="True-0.0", n=200, offset=0.0)
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
    d = @rsubset(effs_overlap, :offset in set)
    u = @rsubset(unique_overlap, :offset in set)

    # Select results for bands
    res = @rsubset(overlap_bands, :offset in set)
    var = :offset

    # Figure options
    f = Figure(; size=(800, 750))
    rev = true

    # Overlap panel
    function make_overlap_ax!(g1, g3; d=d, u=u, l1, rev=rev)
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
            color=:overlap_sign,
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
        ax3 = Axis(g3[1, 1]; xticks=([0.5, 1.5, 2.5], ["Positive", "Overlap", "Negative"]))
        sortedmeasures = reverse(first.(pal))
        m34 = mapping(
            :variable => sorter(sortedcomps),
            [1];
            stack=:countmeasure => sorter(sortedmeasures),
            bar_labels=:label => verbatim,
            color=:countmeasure,
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

        return (g1, g3)
    end

    # Bands panel
    function make_overlap_bands!(f; res=res, var=var, rev=!rev, title="", addlines=false)
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
        f[g4l, 2]; res=@rsubset(res, :offset >= 0.0), title="Overestimation"
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
    Label(f[end:end, end, Bottom()], string(pmax_opt); font=:bold)
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

## Mix overlap & efficiency sign

begin
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    d = @rsubset(effs_overlap, :offset in set)
    u = @rsubset(unique_overlap, :offset in set)

    # Select results for bands
    res_bands = @rsubset(within_bands, :offset >= -0.5, :offset <= 0.5)
    var = :offset

    # Figure options
    f = Figure(; size=(800, 750))
    rev = true

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

##
## Check the distribution of range efficiencies
#=

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

=#
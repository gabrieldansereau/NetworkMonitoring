# using DrWatson
# @quickactivate :NetworkMonitoring

include("include.jl") # see note regarding why we cannot use the module

# Load data
# pmax_opt = "n_at_pmax4"
effs_estimations = CSV.read(datadir("efficiency_estimations.csv"), DataFrame)

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
            r.rmse = max_row.rmse[1]
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

## Plots results

# Visualize
begin
    # Select results for comparison
    set = collect(-0.5:0.1:0.5)
    res_comps = @rsubset(within_comps, :offset in set)
    res_summary = @rsubset(
        unique_comps, :offset in set, :countmeasure in ["positive", "negative", "overlap"]
    )

    # Select results for bands
    # res_overlap = @rsubset(overlap_bands, :offset in collect(-0.5:0.04:0.5))
    res_overlap = @rsubset(overlap_bands, :offset >= -0.5, :offset <= 0.5)
    res_bands = @rsubset(within_bands, :offset >= -0.5, :offset <= 0.5)
    var = :offset

    # Set colour palette
    pal = Dict(
        "negative" => Makie.wong_colors()[2],
        "overlap" => Makie.wong_colors()[4],
        "positive" => Makie.wong_colors()[3],
    )
    scl = scales(; Color=(; palette=[k => v for (k, v) in pal]))
end

# Bands panel
begin
    function make_bands_ax!(
        f; res=res_bands, var=var, rev=false, colour=true, arrowlabels=false
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
            col1 = Makie.wong_colors()[2]
            col2 = Makie.wong_colors()[1]
            lab1 = "Higher efficiency"
            lab2 = "Lower efficiency"
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
        lines!(ax, x, med; label="Median", color=col2)

        # Common options
        vlines!(ax, 0.0; linestyle=:solid, color=:lightgrey, linewidth=2.0)
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
    let
        # Figure options
        f = Figure(; size=(900, 450))

        p1 = f[1:5, 1:3]
        ax = make_bands_ax!(p1; res=res_bands, rev=false, colour=true, arrowlabels=true)
        Legend(f[:, end + 1], ax, "90% Percentile range"; framevisible=false)
        f
    end
end

# Comparison scatter axis
begin
    function make_comps_ax!(ax; res=res_comps, bands=res_bands)
        Random.seed!(42)
        sc = scatter!(
            ax,
            [o + 0.02 * (rand() - 0.5) for o in res.offset],
            res.value;
            color=[pal[v] for v in res.overlap],
            label=[
                uppercasefirst(v) => (; color=pal[v], markersize=8) for v in res.overlap
            ],
            markersize=5,
            alpha=0.9,
        )
        # Set axis limits
        ymin = minimum(bands.low)
        ymax = maximum(bands.upp)
        yabsmax = maximum(abs.([ymin, ymax]))
        yoff = 0.1yabsmax
        ylims!(ax, ymin - yoff, ymax + yoff)
        return sc
    end
    let
        f = Figure(; size=(600, 300))
        ax = Axis(f[1, 1])
        sc = make_comps_ax!(ax; res=res_comps, bands=res_bands)
        f
    end
end

# Overlap lines and bands
begin
    function make_overlap_bands!(
        f; res=res_overlap, var=var, rev=false, title="", addlines=true
    )
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

        # Common options
        if addlines
            vlines!(ax, 0.0; linestyle=:solid, color=:lightgrey, linewidth=2.0)
        end

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

        return ax
    end
    let
        f = Figure(; size=(600, 300))
        make_overlap_bands!(f; res=res_overlap, var=var, rev=false, title="", addlines=true)
        f
    end
end

# Summary count axis
begin
    function make_summary_ax!(gp; res=res_summary)
        ax0 = Axis(
            gp;
            yticks=([0.5, 1.5, 2.5], ["Negative", "Overlap", "Positive"]),
            xaxisposition=:top,
            xticks=sort(res.offset),
        )
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
        draw!(ax0, data(res) * m34 * v3, scl)
        hidexdecorations!(ax0; ticks=false)
        hideydecorations!(ax0; ticklabels=false, ticks=false)
        hidespines!(ax0)
        return ax0
    end
    let
        f = Figure(; size=(700, 200))
        g = GridLayout(f[:, :])
        ax = make_summary_ax!(g[:, :]; res=res_summary)
        f
    end
end

begin
    # Figure options
    f = Figure(; size=(900, 700))

    # GridLayout
    g1 = GridLayout(f[1:5, 1:3])
    g1l = GridLayout(f[1:5, 4])
    g2 = GridLayout(f[(end + 1):(end + 4), 1:3])
    g2l = GridLayout(f[(end - 3):end, 4])

    # Bands
    ax1 = make_bands_ax!(g1[:, :]; res=res_bands, var=var, rev=false, colour=false)
    # Comparison axis
    sc1 = make_comps_ax!(ax1; res=res_comps, bands=res_bands)
    # Overlap bands
    ax2 = make_overlap_bands!(g2[:, :]; res=res_overlap, rev=false, title="")

    # Legends
    leg_opt = (; framevisible=false, merge=true, halign=:left)
    Legend(g1l[:, :], ax1, "Comparison values"; leg_opt...)
    Legend(g2l[:, :], ax2, "Comparison sign"; leg_opt..., tellwidth=false)

    # Align axes
    linkxaxes!(ax1, ax2)

    # Add subpanel labels
    lab_opt = (; halign=:left, font=:bold)
    Label(
        g1[1, 1, Top()],
        "a) Comparison of efficiencies between range estimations";
        lab_opt...,
        padding=(-65, 0, 30, 0),
    )
    Label(
        g2[1, 1, Top()],
        "b) Proportion of simulations per comparison sign";
        lab_opt...,
        padding=(-65, 0, 15, 0),
    )

    # Add over-under info labels
    for g in [g1]
        opt = (; halign=:center, fontsize=12, tellheight=false, tellwidth=false)
        Label(g[1, 1, Top()], "⬅️ Underestimation"; padding=(-200, 0, -15, 0), opt...)
        Label(g[1, 1, Top()], "Overestimation ➡️"; padding=(200, 0, -15, 0), opt...)
    end

    # Optionally add summary panel
    addsummary = false
    if addsummary
        # g0 = GridLayout(f[-1:0, 1:3])
        g0 = GridLayout(f[(end + 1):(end + 2), 1:3])
        ax0 = make_summary_ax!(g0[:, :]; res=res_summary)
        vlines!(ax0, [0.0]; linestyle=:solid, color=:lightgrey, linewidth=2.0)
        linkxaxes!(ax1, ax2, ax0)
    end

    # Figure
    save(plotsdir("ranges_overlap.png"), f)
    f
end

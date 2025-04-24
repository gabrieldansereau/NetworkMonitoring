# using DrWatson
# @quickactivate :NetworkMonitoring

using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFramesMeta
using DrWatson
using Statistics

# Load results
param_grid = CSV.read(datadir("param_grid.csv"), DataFrame)

# Stack all proportion variables
param_stack = stack(
    param_grid,
    filter(startswith("prop"), names(param_grid)),
    value_name=:prop
)

# Reorder variables by expected proportion
order1 = Dict(
    "prop_monitored_sp" => 1,
    "prop_possible_int" => 2,
    "prop_realized_int" => 3,
    "prop_detected_int" => 4,
)
sort!(param_stack, order(:variable, by=x -> order1[x]))

# Rename variables
renamed = Dict(
    "prop_monitored_sp" => "Monitored sp",
    "prop_possible_int" => "Possible int",
    "prop_realized_int" => "Realized int",
    "prop_detected_int" => "Detected int",
)
@rtransform!(param_stack, :variable = renamed[:variable])

## Proportion results

# Common plot
fig = data(param_stack) *
    visual(
        RainClouds;
        markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element",
        layout=:refmethod => renamer("global" => "Global reference", "metawebify" => "Per-element reference"),
    ) |> x ->
    draw(x,
        axis=(; yticks=(0.0:0.25:1.0), xticks=(0:25:100)),
        legend=(; framevisible=false),
        figure=(; size=(700,450))
    )
save(plotsdir("nbon.png"), fig; px_per_unit=2.0)

fig = filter(:refmethod => ==("global"), param_stack) |> x ->
    filter(:variable => in(["Monitored sp"]), x) |> x ->
    # filter(:variable => in(["Monitored sp", "Possible int"]), x) |> x ->
    # filter(:variable => in(["Monitored sp", "Possible int", "Realized int"]), x) |> x ->
    data(x) *
    visual(
        RainClouds;
        markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element",
        # layout=:refmethod => renamer("global" => "Global reference", "metawebify" => "Per-element reference"),
    ) |> x ->
    draw(x,
        axis=(; yticks=(0.0:0.25:1.0), xticks=(0:25:100), limits=((-4, 104), (-0.05, 1.05))),
        legend=(; framevisible=false),
        # figure=(; size=(700,450))
    )
save(plotsdir("nbon_global_1.png"), fig; px_per_unit=2.0)

# Log scale
fig = data(param_stack) *
    visual(
        RainClouds;
        markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => log => "Number of sites in BON (log scale)",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element",
        layout=:refmethod => renamer("global" => "Global reference", "metawebify" => "Per-element reference"),
    ) |> x ->
    draw(x,
        axis=(; yticks=(0.0:0.25:1.0)),
        legend=(; framevisible=false),
        figure=(; size=(700,450))
    )
save(plotsdir("nbon_logx.png"), fig; px_per_unit=2.0)

## Betadiversity results

# Load results
beta_res = CSV.read(datadir("betadiversity.csv"), DataFrame)
filter!(:refmethod => ==("global"), beta_res)

# Summarize
# beta_gp = @chain begin
#     beta_res
#     groupby([:nbon, :nrep])
#     @combine(
#         :βos_possible = median(:βos_possible),
#         :βos_realized = median(:βos_realized),
#         :βos_detected = median(:βos_detected),
#     )
# end
@transform!(beta_res, :monitored_sp = NaN)

# Plot
fig = data(beta_res) *
    visual(
        RainClouds;
        markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        [:monitored_sp, :βos_possible, :βos_realized, :βos_detected] .=> "βos' (monitored vs real metaweb)";
        color=dims(1) =>
            renamer(["Monitored sp", "Possible int", "Realized int", "Detected int"]) =>
            "Sampled Element",
        # color=:variable => presorted => "Sampled element",
        # layout=:refmethod => renamer("global" => "Global reference", "metawebify" => "Per-element reference"),
    ) |> x ->
    draw(x,
        axis=(; yticks=(0.0:0.25:1.0), xticks=(0:25:100)),
        legend=(; framevisible=false),
        # figure=(; size=(700,450))
    )
save(plotsdir("betadiversity.png"), fig; px_per_unit=2.0)

## Median & quantiles

# Load results
beta_meds = CSV.read(datadir("betadiversity_med_quant.csv"), DataFrame)
# filter!(:refmethod => ==("global"), beta_meds)

# Plot
fig = data(beta_meds) *
    visual(
        RainClouds;
        markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        [:βos_upp_possible :βos_med_possible :βos_low_possible;
         :βos_upp_realized :βos_med_realized :βos_low_realized;
         :βos_upp_detected :βos_med_detected :βos_low_detected] .=> "βos' at monitored sites";
        color=dims(2) => renamer(["upp", "med", "low"]) => "Quantile",
        row=dims(1) => renamer(["possible", "realized", "detected"]),
    )
fig_link = draw(fig,
        axis=(; xticks=(0:25:100)),
        legend=(; framevisible=false),
    )
fig_unlink = draw(fig,
        axis=(; xticks=(0:25:100)),
        facet=(; linkyaxes=:none),
        legend=(; framevisible=false),
    )
save(plotsdir("betadiversity_meds_link.png"), fig_link; px_per_unit=2.0)
save(plotsdir("betadiversity_meds_unlink.png"), fig_unlink; px_per_unit=2.0)


# Summarize
beta_gp = @chain begin
    beta_meds
    groupby([:nbon])
    @combine(
        :βos_med_possible = median(:βos_med_possible),
        :βos_med_realized = median(:βos_med_realized),
        :βos_med_detected = median(:βos_med_detected),
    )
    @transform!(:monitored_sp = NaN)
end

# Plot
fig = data(beta_gp) *
    visual(
        RainClouds;
        markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        [:monitored_sp, :βos_med_possible, :βos_med_realized, :βos_med_detected] .=> "Median of βOS' at monitored sites";
        color=dims(1) =>
            renamer(["Monitored sp", "Possible int", "Realized int", "Detected int"]) =>
            "Sampled Element",
        # color=:variable => presorted => "Sampled element",
        # layout=:refmethod => renamer("global" => "Global reference", "metawebify" => "Per-element reference"),
    ) |> x ->
    draw(x,
        axis=(; yticks=(0.0:0.25:1.0), xticks=(0:25:100)),
        legend=(; framevisible=false),
        # figure=(; size=(700,450))
    )
save(plotsdir("betadiversity_meds_meds.png"), fig; px_per_unit=2.0)

## βOS' with monitored metaweb

# Load results
# beta_meds = CSV.read(datadir("betadiversity.csv"), DataFrame)
all_res = wload("./data/all_res.jld2")
beta_mon = all_res["all_res"]
select!(beta_mon, :nbon, r"^βos_mon")
inds_nothing = findall(isnothing, beta_mon.βos_mon_low_detected)
deleteat!(beta_mon, inds_nothing)

# Plot
fig = data(beta_mon) *
    visual(
        RainClouds;
        markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        [:βos_mon_med_possible;
         :βos_mon_med_realized;
         :βos_mon_med_detected] .=> "βos' at monitored sites (vs monitored metaweb)";
        # color=dims(2) => renamer(["upp", "med", "low"]) => "Quantile",
        row=dims(1) => renamer(["possible", "realized", "detected"]),
    ) |> x ->
    draw(x,
        axis=(; xticks=(0:25:100)),
        legend=(; framevisible=false),
    )
save(plotsdir("betadiversity_mon.png"), fig; px_per_unit=2.0)

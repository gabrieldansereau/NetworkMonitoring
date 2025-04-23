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

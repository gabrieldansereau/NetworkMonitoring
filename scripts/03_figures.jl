# using DrWatson
# @quickactivate :NetworkMonitoring

using AlgebraOfGraphics
using CairoMakie
using CSV
using DataFramesMeta
using DrWatson
using Statistics

update_theme!(; CairoMakie=(; px_per_unit=2.0))

## Proportion results

# Load & prep proportion results
function load_and_prep(file)
    # Load results
    param_grid = CSV.read(file, DataFrame)

    # Stack all proportion variables
    param_stack = stack(
        param_grid, filter(startswith("prop"), names(param_grid)); value_name=:prop
    )

    # Reorder variables by expected proportion
    order1 = Dict(
        "prop_monitored_sp" => 1,
        "prop_possible_int" => 2,
        "prop_realized_int" => 3,
        "prop_detected_int" => 4,
    )
    sort!(param_stack, order(:variable; by=x -> order1[x]))

    # Rename variables
    renamed = Dict(
        "prop_monitored_sp" => "Monitored sp",
        "prop_possible_int" => "Possible int",
        "prop_realized_int" => "Realized int",
        "prop_detected_int" => "Detected int",
    )
    @rtransform!(param_stack, :variable = renamed[:variable])
end
param_stack = load_and_prep(datadir("param_grid.csv"))

# Proportion results plot
fig =
    data(param_stack) *
    visual(
        RainClouds; markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element",
        layout=:refmethod => renamer(
            "global" => "Global reference", "metawebify" => "Per-element reference"
        ),
    ) |>
    x -> draw(
        x;
        axis=(; yticks=(0.0:0.25:1.0), xticks=(0:25:100)),
        legend=(; framevisible=false),
        figure=(; size=(700, 450)),
    )
save(plotsdir("nbon.png"), fig)

## Random sampling comparison

# Load
param_rand = load_and_prep(datadir("param_grid-random.csv"))
@transform!(param_rand, :sampling = "random")
@transform!(param_stack, :sampling = "balanced")
append!(param_rand, param_stack)

# Compare with scatter plot
fig =
    data(filter(:refmethod => ==("global"), param_rand)) *
    visual(
        RainClouds; markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        :prop => "Proportion of sampled elements";
        color=:sampling => "Sampling type",
        layout=:variable => presorted,
    ) |>
    x -> draw(
        x;
        axis=(; yticks=(0.0:0.25:1.0), xticks=(0:25:100)),
        legend=(; framevisible=false),
        figure=(; size=(700, 450)),
    )
save(plotsdir("random.png"), fig)

# Compare with violin plot
fig =
    data(param_rand) *
    visual(Violin) *
    mapping(
        :variable => presorted => "Sampled element",
        :prop => "Proportion of sampled elements";
        color=:sampling => "Sampling type",
        dodge=:sampling,
        layout=:refmethod => renamer(
            "global" => "Global reference", "metawebify" => "Per-element reference"
        ),
    ) |>
    x -> draw(
        x;
        axis=(; yticks=(0.0:0.25:1.0), xticklabelrotation=pi / 4),
        legend=(; framevisible=false),
        figure=(; size=(700, 450)),
    )
save(plotsdir("random_violin.png"), fig)

## Betadiversity results

# Load results
beta_res = CSV.read(datadir("betadiversity.csv"), DataFrame)
filter!(:refmethod => ==("global"), beta_res)

# Summarize
@transform!(beta_res, :monitored_sp = NaN)

# Plot
fig =
    data(beta_res) *
    visual(
        RainClouds; markersize=4, jitter_width=0.0, clouds=nothing, plot_boxplots=false
    ) *
    mapping(
        :nbon => "Number of sites in BON",
        [:monitored_sp, :βos_possible, :βos_realized, :βos_detected] .=>
            "βos' (monitored vs real metaweb)";
        color=dims(1) =>
            renamer(["Monitored sp", "Possible int", "Realized int", "Detected int"]) => "Sampled Element",
        # color=:variable => presorted => "Sampled element",
        # layout=:refmethod => renamer("global" => "Global reference", "metawebify" => "Per-element reference"),
    ) |>
    x -> draw(
        x;
        axis=(; yticks=(0.0:0.25:1.0), xticks=(0:25:100)),
        legend=(; framevisible=false),
        # figure=(; size=(700,450))
    )
save(plotsdir("betadiversity.png"), fig; px_per_unit=2.0)

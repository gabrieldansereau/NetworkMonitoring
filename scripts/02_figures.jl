include("00_include.jl")

# Load results
param_grid = CSV.read("data/param_grid.csv", DataFrame)

# Stack all proportion variables
param_stack = stack(
    param_grid,
    filter(startswith("prop"), names(param_grid)),
    value_name=:prop
)

# Reorder variables by expected proportion
order1 = Dict(
    "prop_observed_sp" => 1,
    "prop_possible_int" => 2,
    "prop_realized_int" => 3,
    "prop_detected_int" => 4,
)
sort!(param_stack, order(:variable, by=x -> order1[x]))

# Rename variables
renamed = Dict(
    "prop_observed_sp" => "Observed sp",
    "prop_possible_int" => "Possible int",
    "prop_realized_int" => "Realized int",
    "prop_detected_int" => "Detected int",
)
@rtransform!(param_stack, :variable = renamed[:variable])

## Plot results

# Common jitter plot
jitterplot = data(param_stack) *
    visual(
        RainClouds;
        markersize=8, jitter_width=0.5, clouds=nothing, plot_boxplots=false
    )

# Connectance
fig = jitterplot *
    mapping(
        :C_exp => nonnumeric => "Expected connectance",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element"
    ) |> draw
save("figures/connectance.png", fig; px_per_unit=2.0)

# RA_sigma
fig = jitterplot *
    mapping(
        :ra_sigma => nonnumeric => "Realized abundance sigma value",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element"
    ) |> draw
save("figures/ra_sigma.png", fig; px_per_unit=2.0)

# RA_scaling
fig = jitterplot *
    mapping(
        :ra_scaling => nonnumeric => "Realized abundance scaling value",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element"
    ) |> draw
save("figures/ra_scaling.png", fig; px_per_unit=2.0)

# Energy
fig = jitterplot *
    mapping(
        :energy_NFL => nonnumeric => "Energy - Neutrally forbidden links",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element"
    ) |> draw
save("figures/energy_nfl.png", fig; px_per_unit=2.0)

# Autocorrelation
fig = jitterplot *
    mapping(
        :H_nlm => nonnumeric => "H - Neutral Landscape autocorrelation",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element"
    ) |> draw
save("figures/h_nlm.png", fig; px_per_unit=2.0)

# Number of sampled sites
fig = jitterplot *
    mapping(
        :nbon => nonnumeric => "Number of sites in BON",
        :prop => "Proportion of sampled elements";
        color=:variable => presorted => "Sampled element"
    ) |> draw
save("figures/nbon.png", fig; px_per_unit=2.0)



# Facets?
param_stack2 = stack(param_stack, [:C_exp, :ra_sigma], variable_name=:xvar)
data(param_stack2) *
    mapping(:value => nonnumeric, :prop; color=:variable, layout=:xvar) *
    visual(
        RainClouds;
        markersize=10, jitter_width=0.01, clouds=nothing, plot_boxplots=false
    ) |> x -> draw(x; facet = (; linkxaxes = :none))

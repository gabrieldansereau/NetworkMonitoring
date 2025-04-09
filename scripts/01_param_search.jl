include("00_include.jl")

Random.seed!(42)
SpeciesInteractionSamplers.INTERACTIVE_REPL = false
set_params = false

param_grid = allcombinations(
    DataFrame,
    ns = [100],
    nsites = [100],
    C_exp = [0.1, 0.2, 0.3],
    ra_sigma = [0.6, 1.2, 1.8],
    ra_scaling = [25.0, 50.0, 75.0],
    energy_NFL = [1000, 10_000, 100_000],
    H_nlm = [0.25, 0.5, 0.75],
    nbon = [25, 50, 100],
    prop_observed_sp = missings(Float64, 1),
    prop_detected_int = missings(Float64, 1),
    prop_realized_int = missings(Float64, 1),
    prop_possible_int = missings(Float64, 1),
)

@showprogress for i in 1:nrow(param_grid)
    r = param_grid[i, :]
    global ns = r.ns
    global nsites = r.nsites
    global C_exp = r.C_exp
    global ra_sigma = r.ra_sigma
    global ra_scaling = r.ra_scaling
    global energy_NFL = r.energy_NFL
    global H_nlm = r.H_nlm
    global nbon = r.nbon

    include("main.jl")

    r.prop_observed_sp = round(prop_observed_sp; sigdigits=3)
    r.prop_detected_int = round(prop_detected_int; sigdigits=3)
    r.prop_realized_int = round(prop_realized_int; sigdigits=3)
    r.prop_possible_int = round(prop_possible_int; sigdigits=3)
end
param_grid

CSV.write("param_grid.csv", param_grid)
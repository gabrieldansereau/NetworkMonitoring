using Distributed
addprocs(8)
@everywhere using DrWatson

@everywhere begin
    @quickactivate :NetworkMonitoring

    # Load main script
    include("main.jl")
end

# Define script options
Random.seed!(42)
const SpeciesInteractionSamplers.INTERACTIVE_REPL = false

# Define parameters to explore
const params = Dict(
    :ns => [75],
    :nsites => [100],
    :C_exp => [0.2],
    :ra_sigma => [1.2],
    :ra_scaling => [50.0],
    :energy_NFL => [50_000],
    :H_nlm => [0.5],
    :nbon => collect(1:100),
    :nrep => collect(1:20),
    :refmethod => ["metawebify", "global"]
)
const dicts = dict_list(params)

# Run for all combinations
@showprogress @distributed for d in dicts
    res = main(d)
    d2 = merge!(d, res)
    tagsave(datadir("sim", savename(d, "jld2")), tostringdict(d2))
end

# Test load (with gitcommit)
# wload(datadir("sim", readdir(datadir("sim"))[1]))

# Collect & export results
param_grid = collect_results(datadir("sim"))

# Export results
CSV.write(datadir("param_grid.csv", param_grid))

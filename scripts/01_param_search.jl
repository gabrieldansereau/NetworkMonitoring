using DrWatson
@quickactivate :NetworkMonitoring

# Define script options
Random.seed!(42)
const SpeciesInteractionSamplers.INTERACTIVE_REPL = false

# Load main script
include("main.jl")
main() # precompile

# Define parameters to explore
const params = Dict(
    :ns => [100],
    :nsites => [100],
    :C_exp => [0.1, 0.2, 0.3],
    :ra_sigma => [0.6, 1.2, 1.8],
    :ra_scaling => [25.0, 50.0, 75.0],
    :energy_NFL => [1000, 10_000, 100_000],
    :H_nlm => [0.25, 0.5, 0.75],
    :nbon => [25, 50, 100],
    :refmethod => ["metawebify"]
)
const dicts = dict_list(params)

# Run for all combinations
function run_search!()
    @showprogress for d in dicts
        res = main(d)
        merge!(d, res)
        tagsave(datadir("sim", savename(d, "jld2")), tostringdict(d))
    end
end
run_search!()

# Test load (with gitcommit)
# wload(datadir("sim", readdir(datadir("sim"))[1]))

# Collect & export results
param_grid = collect_results(datadir("sim"))

# Export results
CSV.write(datadir("param_grid.csv", param_grid))

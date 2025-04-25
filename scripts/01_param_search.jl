using Distributed
addprocs(4)
@everywhere using DrWatson

@everywhere begin
    @quickactivate :NetworkMonitoring
end

# Define script options
Random.seed!(42)
const SpeciesInteractionSamplers.INTERACTIVE_REPL = false

## Run parameter search

# Define parameters to explore
const params = Dict(
    :nbon => collect(1:100),
    :nrep => collect(1:20),
    :refmethod => ["metawebify", "global"],
    :res => :monitored,
)
const dicts = dict_list(params)

# Run for all combinations
@showprogress @distributed for d in dicts
    res = runsim(; d...)
    d2 = merge(d, res)
    tagsave(datadir("sim", savename(d, "jld2")), tostringdict(d2))
end

# Test load (with gitcommit)
# wload(datadir("sim", readdir(datadir("sim"))[1]))

## Extract results

# Collect results
param_grid = collect_results(datadir("sim"))
select!(param_grid, Not(:path))

# Sort results
select!(
    param_grid,
    [:nbon, :nrep, :refmethod,
     :prop_monitored_sp, :prop_possible_int, :prop_realized_int, :prop_detected_int]
)
sort!(param_grid, [:refmethod, :nbon, :nrep])

# Remove duplicates?
unique!(param_grid, filter(!startswith("prop"), names(param_grid)))

# Export results
CSV.write(datadir("param_grid.csv"), param_grid)

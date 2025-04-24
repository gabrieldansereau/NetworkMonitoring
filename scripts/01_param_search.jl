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
const _fixed_params = Symbol[]
for (k,v) in params
    if length(v) == 1
        push!(_fixed_params, k)
    end
end

# Run for all combinations
@showprogress @distributed for d in dicts
    res = runsim(d; res=:monitored)
    d2 = merge(d, res)
    tagsave(datadir("sim", savename(d, "jld2"; ignores=_fixed_params)), tostringdict(d2))
end

# Test load (with gitcommit)
# wload(datadir("sim", readdir(datadir("sim"))[1]))

## Extract results

# Collect results
param_grid = collect_results(datadir("sim"))
select!(param_grid, Not([:path, _fixed_params...]))

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

using Distributed

# Instantiate & precompile everywhere
@everywhere begin
    using Pkg; Pkg.activate(@__DIR__)
    Pkg.instantiate(); Pkg.precompile
end

# Load environment
@everywhere begin
    using DrWatson
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
)
const output = :prop # :monitored or :prop
const dicts = dict_list(params)

# Run for all combinations
@showprogress @distributed for d in dicts
    try
        res = runsim(; output=output, d...)
        d2 = merge(d, res)
        tagsave(datadir("sim-$output", savename(d, "jld2")), tostringdict(d2))
        true
    catch e
        false
    end
end

# Test load (with gitcommit)
# wload(readdir(datadir("sim-$output"), join=true)[1]))

## Extract proportion results

if output == :prop

    # Collect results
    param_grid = collect_results(datadir("sim-$output"))
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

end

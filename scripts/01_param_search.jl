# Load environment
using DrWatson
@quickactivate :NetworkMonitoring

# Define script options
Random.seed!(42)
const SpeciesInteractionSamplers.INTERACTIVE_REPL = false

## Run parameter search

# Define parameters to explore
STEP = (OUTDIR == "dev" ? 50 : 10)
const params = Dict(
    :nbon => collect([1, STEP:STEP:500...]), :nrep => collect(1:3), :refmethod => ["global"]
)
const output = :prop # :monitored or :prop
const sampler = BON.BalancedAcceptance # or BON.SimpleRandom
const dicts = dict_list(params)

# Set directory to export results
if !(@isdefined OUTDIR)
    const OUTDIR = "dev/sim-$output" # dev (local), sim-prop (remote)
end
mkpath(datadir(OUTDIR))

# Set progress bar display time from environment variable on cluster
const DT = parse(Float64, get(ENV, "PROGRESS_BARS_DT", "0.1"))

# Run for all combinations
function main()
    @showprogress dt = DT for (i, d) in collect(enumerate(dicts))
        Random.seed!(i)
        try
            res = runsim(; output=output, sampler=sampler, d...)
            d2 = merge(d, res)
            tagsave(datadir(OUTDIR, savename(d, "jld2")), tostringdict(d2); warn=false)
            true
        catch e
            false
        end
    end
end
main()

# Test load (with gitcommit)
# wload(readdir(datadir("sim-$output"), join=true)[1]))

## Extract proportion results

if output == :prop

    # Collect results
    param_grid = collect_results(datadir(OUTDIR))
    select!(param_grid, Not(:path))

    # Sort results
    select!(
        param_grid,
        [
            :nbon,
            :nrep,
            :refmethod,
            :prop_monitored_sp,
            :prop_possible_int,
            :prop_realized_int,
            :prop_detected_int,
        ],
    )
    sort!(param_grid, [:refmethod, :nbon, :nrep])

    # Remove duplicates?
    unique!(param_grid, filter(!startswith("prop"), names(param_grid)))

    # Export results
    CSV.write(datadir("param_grid.csv"), param_grid)
end

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

# Why are there more than expect amount?
_params =  names(param_grid) |> filter(!startswith("prop")) |> filter(!=("path"))
_dups = param_grid[nonunique(param_grid, _params; keep=:noduplicates), :]
replace.(_dups.path, datadir("sim") => "") |>
    x -> replace.(x, ".jld2" => "") #|>
    x -> show(DataFrame(path=x), allrows=true)

_d = dicts[1]
savename(_d)
DrWatson.allaccess(_d)

_fixed_params = Symbol[]
@time for (k,v) in params
    if length(v) == 1
        push!(_fixed_params, k)
    end
end
_fixed_params

savename(_d; ignores=_fixed_params)
savename(_d; accesses=_fixed_params)


_oldpaths = replace.(_dups.path, datadir("sim") => "", "/" => "", ".jld2"=>"")
parse_savename.(_oldpaths)[1]

_oldpaths = replace.(param_grid.path, datadir("sim") => "", "/" => "", ".jld2"=>"")
_expected = savename.(dicts)
intersect(_oldpaths, _expected)
filter(in(_expected), _oldpaths)

_oldpaths[1]

wsave("_tmp.jld2", dicts[1])

# Export results
CSV.write(datadir("param_grid.csv", param_grid))

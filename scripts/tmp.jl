using ProgressMeter
using Base.Threads
using DrWatson
using Distributed

OUTDIR = "efficiency"

includet("tmp-pmax.jl")
@showprogress for id in 1:3
    runid(id)
end

@showprogress @threads for id in 1:200
    runid(id)
end

## Add pmax to result files

include("include.jl")

# Check files
files = filter(contains("monitored_estimations"), readdir(datadir("efficiency")))

# Loop all
for id in 1:200
    @info id

    # Load current data
    idp = lpad(id, 3, "0")
    file = datadir("efficiency", "monitored_estimations-$idp.csv")
    df = CSV.read(file, DataFrame)

    # Get pmax
    pfile = datadir("_xtras", "pmax", "pmax-$idp.csv")
    pdf = CSV.read(pfile, DataFrame)

    # Add pmax
    pmax = Dict(r.offset => r.pmax for r in eachrow(pdf))
    @rtransform!(df, :degmax = ismissing(:monitored) ? missing : pmax[:offset])

    # Order columns
    new = @rselect(
        df, :sim = id, :layer, :offset, :nbon, :rep, :monitored, :deg, :degmax, All()
    )

    # Export
    CSV.write(file, new)
end
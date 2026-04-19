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
# weirdos = [1,4,7,12,19,20,21,22,26,29,34,40,46,51,55,61,62,63,64,66,71,78,81,86,88,89,92,93,96,102,105,106,116,118,121,122,125,128,129,134,145,147,148,154,155,160,171,182,183,192,194,196,198,200]
for id in setdiff(1:200, weirdos)
    # for id in 1:200
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
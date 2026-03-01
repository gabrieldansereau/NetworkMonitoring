# Use job id to vary parameters
files = filter(startswith("monitored_samplers"), readdir(datadir("efficiency")))
ids = sort(parse.(Int, replace.(files, "monitored_samplers-" => "", ".csv" => "")))

# Replace
# id = 2
# ids = 4:100
for id in ids
    @info id
    idp = lpad(id, 2, "0")
    f = "monitored_optimized-$idp.csv"
    old = CSV.read(datadir("efficiency", f), DataFrame)

    new = @select(old, :set = "layers", :sp, :type, :sampler, :layer = :sampler, All())
    @transform!(new, :sampler = "UncertaintySampling")

    f_mod = datadir("focal_array", f)
    if isfile(f_mod)
        mod = CSV.read(f_mod, DataFrame)
        @assert names(new) == names(mod)
        # @assert first(new, 1) == first(mod, 1)
        # @assert first(new, 5) == first(mod, 5)
        @assert first(select(new, Not(1:5)), 5) == first(select(mod, Not(1:5)), 5)
    end

    CSV.write(datadir("efficiency", f), new)
end
using DrWatson
@quickactivate :NetworkMonitoring

# Load parameters
param_grid = CSV.read(datadir("param_grid.csv"), DataFrame)
sim_params = filter(!startswith("prop"), names(param_grid))

# Load results
simfiles = readdir(datadir("sim-monitored"); join=true)
filter!(contains("global"), simfiles)
nrows = length(simfiles)
all_res_df = [DataFrame(
    nbon=Int[],
    nrep=Int[],
    refmethod=String[],
    βos_possible=Float64[],
    βos_realized=Float64[],
    βos_detected=Float64[],
) for i in 1:Threads.nthreads()]
@showprogress Threads.@threads for i in eachindex(simfiles)
    f = simfiles[i]

    # Load simulation results
    d = wload(f)

    # Extract metaweb for the simulation
    m = d["m"]

    # Compute βOS
    βres = Dict{String, Any}()
    for n in ["pos", "realized", "detected"]
        networks = d["networks_$n"]
        n = n == "pos" ? "possible" : n

        # βOS' between monitored and real metaweb
        m_monitored = reduce(union, networks)
        βcomponents = betadiversity(βOS, m, m_monitored)
        βosprime = KGL08(βcomponents)
        βres["βos_$n"] = βosprime
    end

    βres["nbon"] = d["nbon"]
    βres["nrep"] = d["nrep"]
    βres["refmethod"] = d["refmethod"]

    # Export results
    push!(all_res_df[Threads.threadid()], βres; cols=:union)
end
all_res = vcat(all_res_df...)

# Export all βOS results
sort!(all_res, [:refmethod, :nbon, :nrep])
CSV.write(datadir("betadiversity.csv"), all_res)
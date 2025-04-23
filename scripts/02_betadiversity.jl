using DrWatson
@quickactivate :NetworkMonitoring

# Load parameters
param_grid = CSV.read(datadir("param_grid.csv"), DataFrame)
sim_params = filter(!startswith("prop"), names(param_grid))

# Load results
all_res = DataFrame()
simfiles = readdir(datadir("sim"); join=true)
@showprogress for f in simfiles
    # Load simulation results
    d = wload(f)

    # Extract metaweb for the simulation
    m = d["m"]

    # Compute βOS' between monitored and real metaweb
    βres = Dict()
    for n in ["pos", "realized", "detected"]
        networks = d["networks_$n"]
        m_monitored = reduce(union, networks)
        βcomponents = betadiversity(βOS, m, m_monitored)
        βosprime = KGL08(βcomponents)
        n = n == "pos" ? "possible" : n
        βres["βos_$n"] = βosprime
    end

    # Export results
    res = DataFrame(merge(βres, d));
    select!(res, sim_params, r"^β")
    append!(all_res, res)
end

# Export all βOS results
sort!(all_res, [:refmethod, :nbon, :nrep])
CSV.write(datadir("betadiversity.csv"), all_res)
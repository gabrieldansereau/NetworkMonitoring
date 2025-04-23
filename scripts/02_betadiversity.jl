using DrWatson
@quickactivate :NetworkMonitoring

# Load parameters
param_grid = CSV.read(datadir("param_grid.csv"), DataFrame)
sim_params = filter(!startswith("prop"), names(param_grid))

# Load results
simfiles = readdir(datadir("sim"); join=true)
nrows = length(simfiles)
all_res_nt = (
    nbon=Vector{Int}(undef, nrows),
    nrep=Vector{Int}(undef, nrows),
    refmethod=Vector{String}(undef, nrows),
    βos_possible=Vector{Float64}(undef, nrows),
    βos_realized=Vector{Float64}(undef, nrows),
    βos_detected=Vector{Float64}(undef, nrows),
)
@showprogress Threads.@threads for i in eachindex(simfiles)
    f = simfiles[i]

    # Load simulation results
    d = wload(f)

    # Extract metaweb for the simulation
    m = d["m"]

    # Compute βOS' between monitored and real metaweb
    βres = Dict{String, Float64}()
    for n in ["pos", "realized", "detected"]
        networks = d["networks_$n"]
        m_monitored = reduce(union, networks)
        βcomponents = betadiversity(βOS, m, m_monitored)
        βosprime = KGL08(βcomponents)
        n = n == "pos" ? "possible" : n
        βres["βos_$n"] = βosprime
    end

    # Export results
    all_res_nt.nbon[i]=d["nbon"]
    all_res_nt.nrep[i]=d["nrep"]
    all_res_nt.refmethod[i]=d["refmethod"]
    all_res_nt.βos_possible[i]=βres["βos_possible"]
    all_res_nt.βos_realized[i]=βres["βos_realized"]
    all_res_nt.βos_detected[i]=βres["βos_detected"]
end
all_res = DataFrame(all_res_nt)

# Export all βOS results
sort!(all_res, [:refmethod, :nbon, :nrep])
CSV.write(datadir("betadiversity.csv"), all_res)
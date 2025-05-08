using DrWatson
@quickactivate :NetworkMonitoring

Random.seed!(42)
update_theme!(
    CairoMakie=(; px_per_unit=2.0),
)

# Set parameters
d = DefaultParams()

# Generate networks using simulations
nets_dict = generate_networks(d)
@unpack detected, realized, pos, metaweb, ranges = nets_dict

## Test focal monitoring for a single species

# Get species with highest degree for test run
_deg, _sp = findmax(degree(metaweb.metaweb))

# Focal monitoring
function focal_monitoring(sp::Symbol; type::Symbol=:pos)
    net = nets_dict[type]
    deg = degree(net.metaweb, sp)
    monitored_deg = 0
    # nbon = 0

    monitored = DataFrame()
    # while monitored_deg < deg && nbon <= 100
    for i in 1:100
        nbon = i
        bon = generate_bon(; nbon=nbon)

        # _monitored = union(monitor(pos, bon)...)
        monitored_int = monitor(interactions, net, bon; makeunique=true)
        monitored_deg = sum(in.(sp, monitored_int))

        info = (sp = sp, type = type, nbon = nbon, deg = deg, monitored = monitored_deg)
        # @info info
        push!(monitored, info)
    end
    return monitored
end
monitored_sp = focal_monitoring(_sp; type=:pos)

# Visualize result
lines(
    monitored_sp.nbon, monitored_sp.monitored;
    axis=(; xlabel="Sites in BON", ylabel="Monitored interactions")
)

## Repeat focal monitoring by network types

# Run for all types
types = [:pos, :realized, :detected]
monitored_vec = Vector{DataFrame}(undef, 3)
for i in eachindex(types)
    monitored_vec[i] = focal_monitoring(_sp; type=types[i])
end
monitored_types = reduce(vcat, monitored_vec)

# Visualize
begin
    fig = Figure()
    ax = Axis(fig[1,1]; xticks=0:25:100)
    for (i, type) in enumerate(types)
        m = filter(:type => ==(type), monitored_types)
        lines!(m.nbon, m.monitored, label=String(type), color = Makie.wong_colors()[i+1])
        # hlines!(unique(m.deg), linestyle = :dash, color = Makie.wong_colors()[i+1])
    end
    fig[1,2] = Legend(fig, ax, "Network Type", framevisible=false)
    fig
end
save(plotsdir("focal_types.png"), fig)


## Repeat with 4 species with different degrees

# Get species to test
degrees = degree(metaweb.metaweb)
spp = sort(collect(degrees), by=x -> x.second, rev=true)[[1, 25, 50, 75]]

# Repeat focal monitoring per species
monitored_vec = Vector{DataFrame}(undef, 4)
for i in eachindex(spp)
    monitored_vec[i] = focal_monitoring(spp[i].first; type=:pos)
end
monitored_spp = reduce(vcat, monitored_vec)
@rtransform!(monitored_spp, :proportion = :monitored / :deg)

# Visualize result
begin
    fig = Figure(size=(600, 600))
    ax1 = Axis(fig[1,1]; ylabel="Monitored interactions", xticks=0:25:100)
    ax2 = Axis(fig[2,1], ylabel="Proportion monitored", xlabel="Number of sites in BON", xticks=0:25:100)
    for (i, sp) in enumerate(unique(monitored_spp.sp))
        m = filter(:sp => ==(sp), monitored_spp)
        # Monitored int
        lines!(ax1, m.nbon, m.monitored, label=string(sp), color = Makie.wong_colors()[i])
        hlines!(ax1, unique(m.deg), linestyle = :dash, alpha = 0.5, color = Makie.wong_colors()[i])
        # Proportion
        lines!(ax2, m.nbon, m.proportion, label=string(sp), color = Makie.wong_colors()[i])
    end
    fig[:,end+1] = Legend(fig, ax1, "Species", framevisible=false)
    fig
end
save(plotsdir("focal_spp.png"), fig)

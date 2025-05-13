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
nets_dict[:possible] = nets_dict[:pos]
@unpack detected, realized, pos, metaweb, ranges = nets_dict

## Test focal monitoring for a single species

# Get species with highest degree for test run
_deg, _sp = findmax(degree(metaweb.metaweb))

# Focal monitoring
function focal_monitoring(sp::Symbol; type::Symbol=:possible, nbons=1:100, sampler=:BalancedAcceptance)
    net = nets_dict[type]
    if type == :detected
        # fix for detected where the metaweb subfield is not reajusted to detected int
        _mw = metawebify(net)
        _no = parse(Int, replace(string(sp), "node_" => ""))
        deg = sum(_mw[_no, :]) + sum(_mw[:, _no]) + sum(_mw[_no, _no])
    else
        deg = degree(net.metaweb, sp)
    end
    # monitored_deg = 0
    # nbon = 0

    monitored = DataFrame()
    # while monitored_deg < deg && nbon <= 100
    for i in nbons
        nbon = i
        if sampler == :BalancedAcceptance
            bon = generate_bon(; nbon=nbon)
        else
            # Extract species range as layer
            sp_range = SDT.SDMLayer(
                occurrence(ranges)[indexin([_sp], ranges.species)...];
                x=(0.0, d.nsites), y=(0.0, d.nsites)
            )

            # Use WeightedBalancedAcceptance sampling
            bon = BON.sample(sampler(nbon), sp_range)
        end

        # _monitored = union(monitor(pos, bon)...)
        monitored_int = monitor(x -> interactions(render(Binary, x)), net, bon; makeunique=true)
        monitored_deg = sum(in.(sp, monitored_int))

        info = (sp = sp, type = type, nbon = nbon, deg = deg, monitored = monitored_deg)
        # @info info
        push!(monitored, info)
    end
    return monitored
end
monitored_sp = focal_monitoring(_sp; type=:possible, nbons=1:100)
monitored_sp = focal_monitoring(_sp; type=:possible, nbons=1:100, sampler=BON.WeightedBalancedAcceptance)

# Visualize result
begin
    f, ax, l = lines(
        monitored_sp.nbon, monitored_sp.monitored;
        label="possible", color = Makie.wong_colors()[2],
        axis=(; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:20:100)
    )
    hlines!(ax, monitored_sp.deg, linestyle=:dash, alpha = 0.5, color = Makie.wong_colors()[2])
    hlines!(ax, _deg, linestyle=:dash, alpha = 0.5, color=:grey, label="metaweb")
    axislegend(position=:rb)
    f
end


## Repeat focal monitoring by network types

Random.seed!(33)

# Run for all types
types = [:possible, :realized, :detected]
monitored_vec = Vector{DataFrame}(undef, length(types))
@showprogress for i in eachindex(types)
    monitored_vec[i] = focal_monitoring(_sp; type=types[i])
end
monitored_types = reduce(vcat, monitored_vec)

# Visualize
begin
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:20:100)
    for (i, type) in enumerate(types)
        m = filter(:type => ==(type), monitored_types)
        lines!(m.nbon, m.monitored, label=String(type), color = Makie.wong_colors()[i+1])
        # hlines!(unique(m.deg), linestyle = :dash, color = Makie.wong_colors()[i+1])
    end
    hlines!(ax, _deg, linestyle=:dash, alpha = 0.5, color=:grey, label="metaweb")
    axislegend(position=:rc)
    fig
end
save(plotsdir("focal_types.png"), fig)

# Re-run for realized and detected
types = [:realized, :detected]
monitored_vec = Vector{DataFrame}(undef, length(types))
for i in eachindex(types)
    monitored_vec[i] = focal_monitoring(_sp; type=types[i], nbons=[1, 500:500:10_000...])
end
monitored_types = reduce(vcat, monitored_vec)

# Visualize result
begin
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:1000:10_000)
    for (i, type) in enumerate(unique(monitored_types.type))
        m = filter(:type => ==(type), monitored_types)
        lines!(m.nbon, m.monitored, label=string(type), color = Makie.wong_colors()[i+2])
        hlines!(unique(m.deg), linestyle = :dash, color = Makie.wong_colors()[i+2])
    end
    hlines!(ax, _deg, linestyle=:dash, alpha = 0.5, color=:grey, label="metaweb")
    axislegend(position = :rb)
    fig
end
save(plotsdir("focal_types2.png"), fig)

## Repeat with 4 species with different degrees

Random.seed!(101)

# Get species to test
degrees = degree(metaweb.metaweb)
spp = sort(collect(degrees), by=x -> x.second, rev=true)[[1, 25, 50, 70]]

# Repeat focal monitoring per species
monitored_vec = Vector{DataFrame}(undef, 4)
@showprogress for i in eachindex(spp)
    monitored_vec[i] = focal_monitoring(spp[i].first; type=:pos)
end
monitored_spp = reduce(vcat, monitored_vec)
@rtransform!(monitored_spp, :proportion = :monitored / :deg)

# Visualize result
begin
    fig = Figure(size=(600, 600))
    ax1 = Axis(fig[1,1]; ylabel="Monitored interactions", xticks=0:20:100)
    ax2 = Axis(fig[2,1],
        ylabel="Proportion monitored", xlabel="Number of sites in BON",
        xticks=0:20:100, yticks=0:0.2:1.0
    )
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

## Explore variations with different sampler

# Run for all types
begin
    Random.seed!(22)
    samplers = [BON.UncertaintySampling, BON.WeightedBalancedAcceptance, BON.BalancedAcceptance, BON.SimpleRandom]
    monitored_vec = Vector{DataFrame}(undef, length(samplers))
    @showprogress for i in eachindex(samplers)
        monitored_vec[i] = focal_monitoring(_sp; type=:realized, sampler=samplers[i], nbons=1:5:500)
        @rtransform!(monitored_vec[i], :sampler = samplers[i])
    end
    monitored_samplers = reduce(vcat, monitored_vec)
end

# Visualize
labels = ["Uncertainty Sampling", "Weighted Balanced Acceptance", "Balanced Acceptance", "Simple Random"]
begin
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for (i, sampler) in enumerate(samplers)
        m = filter(:sampler => ==(sampler), monitored_samplers)
        lines!(m.nbon, m.monitored, label=labels[i], color = Makie.wong_colors()[i+1])
        # hlines!(unique(m.deg), linestyle = :dash, color = Makie.wong_colors()[i+1])
    end
    hlines!(ax, _deg, linestyle=:dash, alpha = 0.5, color=:grey, label="metaweb")
    Legend(fig[:, end+1], ax)
    fig
end
save(plotsdir("focal_samplers.png"), fig)

# Run with replicates
begin
    Random.seed!(22)
    nrep = 5
    samplers = [BON.UncertaintySampling, BON.WeightedBalancedAcceptance, BON.BalancedAcceptance, BON.SimpleRandom]
    monitored_mat = Matrix{DataFrame}(undef, length(samplers), nrep)
    @showprogress for i in eachindex(samplers), j in 1:nrep
        monitored_mat[i,j] = focal_monitoring(_sp; type=:realized, sampler=samplers[i], nbons=1:5:500)
        @rtransform!(monitored_mat[i,j], :sampler = samplers[i])
    end
    monitored_samplers = reduce(vcat, vec(monitored_mat))
end

# Interval bands
bands = @chain begin
    monitored_samplers
    groupby([:sampler, :nbon])
    @combine(
        :low = minimum(:monitored),
        :med = median(:monitored),
        :upp = maximum(:monitored),
    )
end
begin
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for (i, s) in enumerate(unique(bands.sampler))
        b = filter(:sampler => ==(s), bands)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=labels[i], color = Makie.wong_colors()[i+1])
        lines!(b.nbon, b.med, label=labels[i], color = Makie.wong_colors()[i+1])
    end
    # hlines!(ax, _deg, linestyle=:dash, alpha = 0.5, color=:grey, label="metaweb")
    # Legend(fig[1, end+1], ax)
    axislegend(position=:lt, merge=true)
    fig
end
save(plotsdir("focal_samplers_bands.png"), fig)

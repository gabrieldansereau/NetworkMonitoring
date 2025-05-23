using DrWatson
@quickactivate :NetworkMonitoring

using BiodiversityObservationNetworks:
    UncertaintySampling, BalancedAcceptance, WeightedBalancedAcceptance, SimpleRandom

Random.seed!(42)
update_theme!(; CairoMakie=(; px_per_unit=2.0))

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
function focal_monitoring(
    sp::Symbol,
    layer=nothing;
    nrep::Int=1,
    type::Symbol=:possible,
    nbons=1:100,
    sampler=BON.BalancedAcceptance,
    combined=false,
)
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
    @showprogress for rep in 1:nrep
        for i in nbons
            nbon = i
            if isnothing(layer)
                @unpack H_nlm, nsites = d
                layer = SDT.SDMLayer(
                    MidpointDisplacement(H_nlm),
                    (nsites, nsites);
                    x=(0.0, nsites),
                    y=(0.0, nsites),
                )
            end

            # Define monitoring sites on layer given sampler
            bon = BON.sample(sampler(nbon), layer)

            # _monitored = union(monitor(pos, bon)...)
            monitored_int = monitor(
                x -> interactions(render(Binary, x)), net, bon; makeunique=true
            )
            monitored_deg = sum(in.(sp, monitored_int))

            info = (sp=sp, type=type, nbon=nbon, rep=rep, deg=deg, monitored=monitored_deg)
            # @info info
            push!(monitored, info)
        end
    end

    if combined
        monitored = @chain monitored begin
            groupby([:sp, :type, :nbon])
            @combine(
                :low = minimum(:monitored),
                :med = median(:monitored),
                :upp = maximum(:monitored),
                :deg = maximum(:deg)
            )
            rename(:type => :var)
        end
    end

    return monitored
end
monitored_sp = focal_monitoring(_sp; type=:possible, nbons=1:100)
monitored_sp = focal_monitoring(
    _sp; type=:possible, nbons=1:100, sampler=BON.WeightedBalancedAcceptance
)

# Visualize result
begin
    f, ax, l = lines(
        monitored_sp.nbon,
        monitored_sp.monitored;
        label="possible",
        color=Makie.wong_colors()[2],
        axis=(; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:20:100),
    )
    hlines!(ax, monitored_sp.deg; linestyle=:dash, alpha=0.5, color=Makie.wong_colors()[2])
    hlines!(ax, _deg; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb)
    f
end

## Repeat focal monitoring by network types

Random.seed!(33)

# Run for all types
types = [:possible, :realized, :detected]
monitored_vec = Vector{DataFrame}(undef, length(types))
for i in eachindex(types)
    @info "$i: $(types[i])"
    monitored_vec[i] = focal_monitoring(_sp; type=types[i], nrep=5, combined=true)
end
monitored_types = reduce(vcat, monitored_vec)

# Interval bands
labs = Dict{Any,String}(
    :possible => "possible", :realized => "realized", :detected => "detected"
)
cols = Dict(
    :possible => Makie.wong_colors()[2],
    :realized => Makie.wong_colors()[3],
    :detected => Makie.wong_colors()[4],
)
begin
    res = monitored_types
    vars = types
    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:25:100
    )
    for v in vars
        b = filter(:var => ==(v), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=labs[v], color=cols[v])
        lines!(b.nbon, b.med; label=labs[v], color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rc, merge=true)
    fig
end
save(plotsdir("focal_types.png"), fig)

# Re-run for realized and detected
types2 = [:realized, :detected]
monitored_types2 = filter(:var => !(in(types2)), monitored_types)
monitored_vec = Vector{DataFrame}(undef, length(types2))
for i in eachindex(types2)
    @info "$i: $(types[i])"
    monitored_vec[i] = focal_monitoring(
        _sp; type=types2[i], nbons=[1, 500:500:10_000...], nrep=5, combined=true
    )
end
append!(monitored_types2, reduce(vcat, monitored_vec))

# Visualize result
begin
    res = monitored_types2
    vars = types2
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:2000:10_000,
    )
    for v in vars
        b = filter(:var => ==(v), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=labs[v], color=cols[v])
        lines!(b.nbon, b.med; label=labs[v], color=cols[v])
        hlines!(unique(b.deg); linestyle=:dash, color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb, merge=true)
    fig
end
save(plotsdir("focal_types2.png"), fig)

## Repeat with 4 species with different degrees

Random.seed!(101)

# Get species to test
degrees = degree(metaweb.metaweb)
spp = sort(collect(degrees); by=x -> x.second, rev=true)[[1, 25, 50, 70]]

# Repeat focal monitoring per species
monitored_vec = Vector{DataFrame}(undef, 4)
@showprogress for i in eachindex(spp)
    monitored_vec[i] = focal_monitoring(spp[i].first; type=:possible, nrep=5, combined=true)
end
monitored_spp = reduce(vcat, monitored_vec)

# Visualize result
begin
    res = monitored_spp
    fig = Figure(; size=(600, 600))
    ax1 = Axis(fig[1, 1]; ylabel="Monitored interactions", xticks=0:20:100)
    ax2 = Axis(
        fig[2, 1];
        ylabel="Proportion monitored",
        xlabel="Number of sites in BON",
        xticks=0:20:100,
        yticks=0:0.2:1.0,
    )
    for (i, sp) in enumerate(spp)
        sp = sp[1]
        b = filter(:sp => ==(sp), res)
        # Monitored int
        band!(
            ax1,
            b.nbon,
            b.low,
            b.upp;
            alpha=0.4,
            label=string(sp),
            colormap=:seaborn_colorblind,
        )
        lines!(ax1, b.nbon, b.med; label=string(sp), colormap=:seaborn_colorblind)
        hlines!(
            ax1, unique(b.deg); linestyle=:dash, alpha=0.5, colormap=:seaborn_colorblind
        )
        # Proportion
        band!(
            ax2,
            b.nbon,
            b.low ./ b.deg,
            b.upp ./ b.deg;
            alpha=0.4,
            label=string(sp),
            colormap=:seaborn_colorblind,
        )
        lines!(ax2, b.nbon, b.med ./ b.deg; label=string(sp), colormap=:seaborn_colorblind)
    end
    fig[:, end + 1] = Legend(fig, ax1, "Species"; framevisible=false, merge=true)
    fig
end
save(plotsdir("focal_spp.png"), fig)

## Explore variations with different sampler

# Extract species range
_sp_range = SDT.SDMLayer(
    occurrence(ranges)[indexin([_sp], ranges.species)...];
    x=(0.0, d.nsites),
    y=(0.0, d.nsites),
)

# Run with replicates
begin
    Random.seed!(22)
    nrep = 5
    samplers = [UncertaintySampling, WeightedBalancedAcceptance, SimpleRandom]
    monitored_mat = Matrix{DataFrame}(undef, length(samplers), nrep)
    @showprogress for i in eachindex(samplers), j in 1:nrep
        monitored_mat[i, j] = focal_monitoring(
            _sp, _sp_range; type=:realized, sampler=samplers[i], nbons=1:5:500
        )
        @rtransform!(monitored_mat[i, j], :sampler = samplers[i])
    end
    monitored_samplers = reduce(vcat, vec(monitored_mat))
end

# Combine for interval bands
bands = @chain begin
    monitored_samplers
    groupby([:sampler, :nbon])
    @combine(
        :low = minimum(:monitored), :med = median(:monitored), :upp = maximum(:monitored),
    )
end

# Define labels and colors
if !(@isdefined labs)
    labs = Dict{Any,String}()
end
labs[UncertaintySampling] = "Uncertainty Sampling"
labs[WeightedBalancedAcceptance] = "Weighted Balanced Acceptance"
labs[SimpleRandom] = "Simple Random"
if !(@isdefined cols)
    cols = Dict{Any,Any}()
end
cols[UncertaintySampling] = Makie.wong_colors()[2]
cols[WeightedBalancedAcceptance] = Makie.wong_colors()[3]
cols[SimpleRandom] = Makie.wong_colors()[4]

# Plot
begin
    res = bands
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=labs[s], color=cols[s])
        lines!(b.nbon, b.med; label=labs[s], color=cols[s])
    end
    axislegend(; position=:lt, merge=true)
    fig
end
save(plotsdir("focal_samplers_bands.png"), fig)

## Richness-focused sampling & suspicion monitoring

# Extract richness layers
richness_spp = SDT.SDMLayer(sum(occurrence(ranges)); x=(0.0, d.nsites), y=(0.0, d.nsites))
richness_int = SDT.SDMLayer(
    extract(SIN.links, realized); x=(0.0, d.nsites), y=(0.0, d.nsites)
)

# Optimize with UncertaintySampling
omnilabels = ["Richness focused", "Interaction focused"]
begin
    Random.seed!(22)
    nrep = 5
    omnis = [richness_spp, richness_int]
    monitored_mat = Matrix{DataFrame}(undef, length(omnis), nrep)
    @showprogress for (i, omni) in enumerate(omnis), j in 1:nrep
        monitored_mat[i, j] = focal_monitoring(
            _sp, omni; type=:realized, sampler=UncertaintySampling, nbons=1:5:500
        )
        @rtransform!(monitored_mat[i, j], :sampler = omnilabels[i])
    end
    monitored_omni = reduce(vcat, vec(monitored_mat))
end

# Combine with Uncertainty Sampling on species layer
monitored_samplers2 = filter(
    :sampler => in(string.([UncertaintySampling])), monitored_samplers
)
append!(monitored_samplers2, monitored_omni)

# Visualize
labels = ["Uncertainty Sampling", omnilabels...]
bands = @chain begin
    monitored_samplers2
    groupby([:sampler, :nbon])
    @combine(
        :low = minimum(:monitored), :med = median(:monitored), :upp = maximum(:monitored),
    )
end
begin
    res = bands
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for (i, s) in enumerate(unique(res.sampler))
        b = filter(:sampler => ==(s), res)
        band!(
            b.nbon,
            b.low,
            b.upp;
            alpha=0.4,
            label=labels[i],
            color=Makie.wong_colors()[i + 1],
        )
        lines!(b.nbon, b.med; label=labels[i], color=Makie.wong_colors()[i + 1])
    end
    # hlines!(ax, _deg, linestyle=:dash, alpha = 0.5, color=:grey, label="metaweb")
    # Legend(fig[1, end+1], ax)
    axislegend(; position=:lt, merge=true)
    fig
end
save(plotsdir("focal_samplers_richness.png"), fig)
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
deg, sp = findmax(degree(metaweb.metaweb))

# Focal monitoring
function focal_monitoring(
    nets_dict,
    sp::Symbol,
    layer::Union{Nothing,SDT.SDMLayer}=nothing;
    nrep::Int=1,
    nbons::AbstractRange{Int64}=1:100,
    type::Vector{Symbol}=[:possible],
    sampler::Vector{UnionAll}=[BON.BalancedAcceptance],
    combined=false,
    H_nlm=d.H_nlm,
    nsites=d.nsites,
)
    # Define options for simulations
    options_dict = Dict(
        :nrep => collect(1:nrep),
        :nbon => collect(nbons),
        :type => type,
        :sampler => sampler,
    )
    options_list = dict_list(options_dict)

    # Run all options
    monitored = DataFrame()
    @showprogress for opt in options_list
        @unpack nrep, nbon, type, sampler = opt

        # Select network object
        net = nets_dict[type]

        # Get degree of focal species
        if type == :detected
            # fix for detected where the metaweb subfield is not reajusted to detected int
            _mw = metawebify(net)
            _no = parse(Int, replace(string(sp), "node_" => ""))
            deg = sum(_mw[_no, :]) + sum(_mw[:, _no]) + sum(_mw[_no, _no])
        else
            deg = degree(net.metaweb, sp)
        end

        # Create neutral layer for optimization if nore was provided
        if isnothing(layer)
            layer = SDT.SDMLayer(
                MidpointDisplacement(H_nlm),
                (nsites, nsites);
                x=(0.0, nsites),
                y=(0.0, nsites),
            )
        end

        # Generate monitoring sites on layer given sampler
        bon = BON.sample(sampler(nbon), layer)

        # Compute degree for focal species across monitored sites
        monitored_int = monitor(
            x -> interactions(render(Binary, x)), net, bon; makeunique=true
        )
        monitored_deg = sum(in.(sp, monitored_int))

        # Export results
        info = (
            sp=sp,
            type=type,
            sampler=sampler,
            nbon=nbon,
            rep=nrep,
            deg=deg,
            monitored=monitored_deg,
        )
        push!(monitored, info)
    end

    # Combine results across replicates
    if combined
        monitored = @chain monitored begin
            groupby([:sp, :type, :sampler, :nbon])
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
function focal_monitoring(nets_dict, spp::Vector{Symbol}; kw...)
    monitored_vec = Vector{DataFrame}(undef, length(spp))
    for (i, sp) in enumerate(spp)
        @info "Monitoring $sp ($i/$(length(spp))"
        monitored_vec[i] = focal_monitoring(nets_dict, sp; kw...)
    end
    monitored = reduce(vcat, monitored_vec)
    return monitored
end
function focal_monitoring(nets_dict, sp::Symbol, layers::Vector{SDT.SDMLayer}; kw...)
    monitored_vec = Vector{DataFrame}(undef, length(layers))
    for (i, layer) in enumerate(layers)
        @info "Monitoring layer $i/$(length(layers))"
        monitored_vec[i] = focal_monitoring(nets_dict, sp, layer; kw...)
        @rtransform!(monitored_vec[i], :layer = i)
    end
    monitored = reduce(vcat, monitored_vec)
    return monitored
end

# Test run
monitored_sp = focal_monitoring(nets_dict, sp; type=[:possible], nbons=1:100)

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
    hlines!(ax, deg; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb)
    f
end

## Repeat focal monitoring by network types

Random.seed!(33)

# Run for all types
types = [:possible, :realized, :detected]
monitored_types = focal_monitoring(nets_dict, sp; type=types, nrep=5, combined=true)

# Interval bands
labs = Dict{Any,String}(
    :possible => "possible", :realized => "realized", :detected => "detected"
)
cols = Dict{Any,Any}(
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
monitored_types2 = focal_monitoring(
    nets_dict, sp; type=types2, nbons=1:500:10_001, nrep=5, combined=true
)

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
spp = [sp.first for sp in spp]

# Repeat focal monitoring per species
monitored_spp = focal_monitoring(nets_dict, spp; type=[:possible], nrep=5, combined=true)

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
        b = filter(:sp => ==(sp), res)
        # Monitored int
        band!(ax1, b.nbon, b.low, b.upp; alpha=0.4, label=string(sp))
        lines!(ax1, b.nbon, b.med; label=string(sp))
        hlines!(ax1, unique(b.deg); linestyle=:dash, alpha=0.5)
        # Proportion
        band!(ax2, b.nbon, b.low ./ b.deg, b.upp ./ b.deg; alpha=0.4, label=string(sp))
        lines!(ax2, b.nbon, b.med ./ b.deg; label=string(sp))
    end
    fig[:, end + 1] = Legend(fig, ax1, "Species"; framevisible=false, merge=true)
    fig
end
save(plotsdir("focal_spp.png"), fig)

## Explore variations with different sampler

# Extract species range
sp_range = SDT.SDMLayer(
    occurrence(ranges)[indexin([sp], ranges.species)...];
    x=(0.0, d.nsites),
    y=(0.0, d.nsites),
)

# Run with replicates
Random.seed!(22)
samplers = [UncertaintySampling, WeightedBalancedAcceptance, SimpleRandom]
monitored_samplers = focal_monitoring(
    nets_dict,
    sp,
    sp_range;
    type=[:realized],
    sampler=samplers,
    nbons=1:5:500,
    nrep=5,
    combined=true,
)

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
cols[SimpleRandom] = Makie.wong_colors()[1];
cols[UncertaintySampling] = Makie.wong_colors()[2];
cols[WeightedBalancedAcceptance] = Makie.wong_colors()[3];

# Plot
begin
    res = monitored_samplers
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=labs[s], color=cols[s])
        lines!(b.nbon, b.med; label=labs[s], color=cols[s])
    end
    hlines!(ax, deg; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:lt, merge=true)
    fig
end
save(plotsdir("focal_samplers.png"), fig)

## Richness-focused sampling

# Extract richness layers
richness_spp = SDT.SDMLayer(sum(occurrence(ranges)); x=(0.0, d.nsites), y=(0.0, d.nsites))
richness_int = SDT.SDMLayer(
    extract(SIN.links, realized); x=(0.0, d.nsites), y=(0.0, d.nsites)
)

# Extract richness of interacting species for focal species
degree_possible = SDT.SDMLayer(
    extract(x -> degree(render(Binary, x), sp), pos); x=(0.0, d.nsites), y=(0.0, d.nsites)
)
degree_realized = SDT.SDMLayer(
    extract(x -> degree(render(Binary, x), sp), realized);
    x=(0.0, d.nsites),
    y=(0.0, d.nsites),
)

# Optimize with UncertaintySampling
optim = [richness_spp, degree_realized]
optimlabels = ["Species richness", "Realized interactions"]
monitored_optimized = focal_monitoring(
    nets_dict,
    sp,
    optim;
    type=[:realized],
    sampler=[UncertaintySampling],
    nbons=1:5:500,
    nrep=5,
    combined=true,
)
@rtransform!(monitored_optimized, :sampler = optimlabels[:layer])
select!(monitored_optimized, Not(:layer))

# Combine with Uncertainty Sampling on species layer
append!(
    monitored_optimized,
    filter(:sampler => ==(UncertaintySampling), monitored_samplers);
    promote=true,
    cols=:subset,
)

# Set labels and colors
labs[UncertaintySampling] = "Focal species range"
for (i, o) in enumerate(optimlabels)
    labs[o] = o
    cols[o] = Makie.wong_colors()[i + 3]
end

# Plot
begin
    res = monitored_optimized
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=labs[s], color=cols[s])
        lines!(b.nbon, b.med; label=labs[s], color=cols[s])
    end
    hlines!(ax, deg; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend("Optimized on"; position=:lt, merge=true)
    fig
end
save(plotsdir("focal_optimized.png"), fig)

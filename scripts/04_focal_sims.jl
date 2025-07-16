using DrWatson
@quickactivate :NetworkMonitoring

using BiodiversityObservationNetworks:
    UncertaintySampling, BalancedAcceptance, WeightedBalancedAcceptance, SimpleRandom

# Set default parameters for network simulations
d = DefaultParams()

# Generate networks using simulations
Random.seed!(42)
nets_dict = generate_networks(d)
nets_dict[:possible] = nets_dict[:pos]
@unpack detected, realized, pos, metaweb, ranges = nets_dict

## Test focal monitoring for a single species

# Get species with highest degree for test run
deg, sp = findmax(degree(metaweb.metaweb))

# Set number of replicates for short interactive run
if !(@isdefined NREP)
    const NREP = 5
end
# For longer, non-interactive runs, set the number of replicates with:
# julia --project -e 'const NREP = 100; include("scripts/04_focal_sims.jl")'

# Set progress bar display time from environment variable on cluster
const DT = parse(Float64, get(ENV, "PROGRESS_BARS_DT", "0.1"))

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
    @showprogress dt = DT for opt in options_list
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

## Repeat focal monitoring by network types

# Random.seed!(33)

# # Run for all types
# types = [:possible, :realized, :detected]
# monitored_types = focal_monitoring(nets_dict, sp; type=types, nrep=NREP, combined=true)

# # Re-run for realized and detected with more sites in BON
# types2 = [:realized, :detected]
# monitored_types2 = focal_monitoring(
#     nets_dict, sp; type=types2, nbons=1:500:10_001, nrep=NREP, combined=true
# )

# # Export
# CSV.write(datadir("monitored_types.csv"), monitored_types)
# CSV.write(datadir("monitored_types2.csv"), monitored_types2)

# ## Repeat with 4 species with different degrees

Random.seed!(101)

# Get species to test
degrees = degree(metaweb.metaweb)
spp = sort(collect(degrees); by=x -> x.second, rev=true)[[1, 25, 50, 70]]
spp = [sp.first for sp in spp]

# Repeat focal monitoring per species
monitored_spp = focal_monitoring(
    nets_dict, spp; type=[:realized], nrep=NREP, nbons=1:5:500, combined=true
)

# Export
CSV.write(datadir("monitored_spp.csv"), monitored_spp)

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
    nrep=NREP,
    combined=true,
)

# Export
CSV.write(datadir("monitored_samplers.csv"), monitored_samplers)
SDT.SimpleSDMLayers.save(datadir("layer_sp_range.tiff"), sp_range)

## Richness-focused sampling

# Extract richness layers
richness_spp = SDT.SDMLayer(sum(occurrence(ranges)); x=(0.0, d.nsites), y=(0.0, d.nsites))
richness_int = SDT.SDMLayer(
    extract(SIN.links, realized); x=(0.0, d.nsites), y=(0.0, d.nsites)
)
richness_pos = SDT.SDMLayer(extract(SIN.links, pos); x=(0.0, d.nsites), y=(0.0, d.nsites))

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
    nrep=NREP,
    combined=true,
)
@rtransform!(monitored_optimized, :sampler = optimlabels[:layer])

# Combine with Uncertainty Sampling on focal species layer
@chain begin
    monitored_samplers
    filter(:sampler => ==(UncertaintySampling), _)
    @transform!(:sampler = "Focal species range")
    append!(monitored_optimized, _; promote=true, cols=:subset)
end

# Export
CSV.write(datadir("monitored_optimized.csv"), monitored_optimized)
SDT.SimpleSDMLayers.save(datadir("layer_richness_spp.tiff"), richness_spp)
SDT.SimpleSDMLayers.save(datadir("layer_richness_int.tiff"), richness_int)
SDT.SimpleSDMLayers.save(datadir("layer_richness_pos.tiff"), richness_pos)
SDT.SimpleSDMLayers.save(datadir("layer_degree_realized.tiff"), degree_realized)
SDT.SimpleSDMLayers.save(datadir("layer_degree_possible.tiff"), degree_possible)

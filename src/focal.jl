function generate_focal_simulation(d; seed=rand(1:10_000))
    # Generate networks using simulations
    Random.seed!(seed)
    nets_dict = generate_networks(d)
    nets_dict[:possible] = nets_dict[:pos]
    @unpack detected, realized, pos, metaweb, ranges = nets_dict

    # Generate probabilistic ranges
    # We need to run operations in the same order
    Random.seed!(seed)
    _ = generate(SIS.NicheModel(d.ns, d.C_exp)) # only to match random seed in next function
    probranges = generate(
        AutocorrelatedProbabilisticRange(; dims=(d.nsites, d.nsites)), d.ns
    )

    # Get species with highest degree
    deg, sp = findmax(degree(metaweb.metaweb))

    # Extract species range
    sp_range = SDT.SDMLayer(
        occurrence(ranges)[indexin([sp], ranges.species)...];
        x=(0.0, d.nsites),
        y=(0.0, d.nsites),
    )

    # Create species range mask
    sp_mask = SDT.nodata(sp_range, 0)

    # Extract richness layers
    richness_spp = SDT.SDMLayer(
        sum(occurrence(ranges)); x=(0.0, d.nsites), y=(0.0, d.nsites)
    )
    richness_int = SDT.SDMLayer(
        extract(SIN.links, realized); x=(0.0, d.nsites), y=(0.0, d.nsites)
    )
    richness_pos = SDT.SDMLayer(
        extract(SIN.links, pos); x=(0.0, d.nsites), y=(0.0, d.nsites)
    )

    # Extract richness of interacting species for focal species
    degree_possible = SDT.SDMLayer(
        extract(x -> degree(render(Binary, x), sp), pos);
        x=(0.0, d.nsites),
        y=(0.0, d.nsites),
    )
    degree_realized = SDT.SDMLayer(
        extract(x -> degree(render(Binary, x), sp), realized);
        x=(0.0, d.nsites),
        y=(0.0, d.nsites),
    )

    # Extract probabilistic range
    probsp_range = SDT.SDMLayer(
        occurrence(probranges)[indexin([sp], probranges.species)...];
        x=(0.0, d.nsites),
        y=(0.0, d.nsites),
    )

    # Extract thresholds
    begin
        Random.seed!(seed)
        thresholds = []
        _ = generate(SIS.NicheModel(d.ns, d.C_exp)) # only to match random seed in next function
        for n in 1:(d.ns)
            ar = AutocorrelatedProbabilisticRange(; dims=(d.nsites, d.nsites))
            H, sz, thres, bin = ar.autocorrelation, ar.dims, rand(ar.threshold), ar.binary
            range_mat = rand(DiamondSquare(H), sz) # kept to match random seed
            push!(thresholds, thres)
        end
        thresholds
    end

    sims_dict = @dict probranges deg sp sp_range sp_mask richness_spp richness_int richness_pos degree_possible degree_realized probsp_range thresholds

    return (nets_dict, sims_dict)
end

function focal_monitoring(
    nets_dict,
    sp::Symbol,
    layer::Union{Nothing,SDT.SDMLayer}=nothing;
    name::String="unnamed",
    nrep::Int=1,
    nbons::AbstractRange{Int64}=1:100,
    type::Vector{Symbol}=[:possible],
    sampler::Vector{UnionAll}=[BON.BalancedAcceptance],
    combined=false,
    H_nlm=DefaultParams().H_nlm,
    nsites=DefaultParams().nsites,
)
    # Define options for simulations
    options_dict = Dict(
        :nrep => collect(1:nrep),
        :nbon => collect(nbons),
        :type => type,
        :sampler => sampler,
    )
    options_list = dict_list(options_dict)

    # Set progress bar display time from environment variable on cluster
    DT = parse(Float64, get(ENV, "PROGRESS_BARS_DT", "0.1"))

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

        # Create neutral layer for optimization if none was provided
        if isnothing(layer)
            layer = SDT.SDMLayer(
                MidpointDisplacement(H_nlm),
                (nsites, nsites);
                x=(0.0, nsites),
                y=(0.0, nsites),
            )
        end

        # Prevent issues with layers without any presence
        if iszero(length(layer))
            monitored_deg = 0
        else
            # Make sure the sampler will work
            len = length(layer)
            if len < nbon && sampler == BON.SimpleRandom
                error("SimpleRandom cannot sample fewer sites than present in layer")
            end
            if len != length(layer.grid) && sampler == BON.WeightedBalancedAcceptance
                error(
                    "WeightedBalancedAcceptance does not work when some pixels have no data"
                )
            end

            # Generate monitoring sites on layer given sampler
            bon = BON.sample(sampler(nbon), layer)

            # Compute degree for focal species across monitored sites
            monitored_int = monitor(
                x -> interactions(render(Binary, x)), net, bon; makeunique=true
            )
            monitored_deg = sum(in.(sp, monitored_int))
        end

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

    # Add simulation set name
    @select!(monitored, :set = name, All())

    return monitored
end
function focal_monitoring(nets_dict, spp::Vector{Symbol}; kw...)
    monitored_vec = Vector{DataFrame}(undef, length(spp))
    for (i, sp) in enumerate(spp)
        @info "Monitoring $sp ($i/$(length(spp)))"
        monitored_vec[i] = focal_monitoring(nets_dict, sp; kw...)
    end
    monitored = reduce(vcat, monitored_vec)
    return monitored
end
function focal_monitoring(
    nets_dict, sp::Symbol, layers::Vector{T}; kw...
) where {T<:SDT.SDMLayer}
    monitored_vec = Vector{DataFrame}(undef, length(layers))
    for (i, layer) in enumerate(layers)
        @info "Monitoring layer $i/$(length(layers))"
        monitored_vec[i] = focal_monitoring(nets_dict, sp, layer; kw...)
        @rtransform!(monitored_vec[i], :layer = i)
    end
    monitored = reduce(vcat, monitored_vec)
    return monitored
end

function summarize_focal(df; id=0)
    if "layer" in names(df)
        cols = [:set, :sp, :type, :sampler, :layer, :nbon]
    else
        cols = [:set, :sp, :type, :sampler, :nbon]
    end
    monitored = @chain df begin
        groupby(cols)
        @combine(
            :low = quantile(:monitored, 0.05),
            :med = median(:monitored),
            :upp = quantile(:monitored, 0.95),
            :deg = maximum(:deg)
        )
        @rtransform(:low = :low / :deg, :med = :med / :deg, :upp = :upp / :deg,)
        @select(:sim = id, All())
    end
    monitored.sampler =
        replace.(
            monitored.sampler,
            "UncertaintySampling" => "Uncertainty Sampling",
            "WeightedBalancedAcceptance" => "Weighted Balanced Acceptance",
            "BalancedAcceptance" => "Balanced Acceptance",
            "SimpleRandomMask" => "Simple Random Mask",
            "SimpleRandom" => "Simple Random",
        )
    return monitored
end

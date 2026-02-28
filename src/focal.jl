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
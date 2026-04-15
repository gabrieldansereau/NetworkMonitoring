using DrWatson
@quickactivate :NetworkMonitoring

using BiodiversityObservationNetworks:
    UncertaintySampling, BalancedAcceptance, WeightedBalancedAcceptance, SimpleRandom

# Set directory to export results
if !(@isdefined OUTDIR)
    const OUTDIR = "dev" # dev (local), focal_array or efficiency (remote)
end
mkpath(datadir(OUTDIR))

function runid(id)
    # Use job id to random results
    # id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    idp = lpad(id, 3, "0")

    # Generate networks using simulations
    d = DefaultParams()
    seed = id * 42
    nets_dict, sims_dict = generate_focal_simulation(d; seed=seed)
    @unpack detected, realized, pos, metaweb, ranges = nets_dict
    @unpack probranges, deg, sp, sp_range, sp_mask = sims_dict
    @unpack richness_spp, richness_int, richness_pos, degree_possible, degree_realized =
        sims_dict
    @unpack probsp_range, thresholds = sims_dict

    ## Range estimation

    # Extract the current threshold
    threshold = thresholds[indexin([sp], probranges.species)...]

    # Define the estimation errors to explore (percentage of range)
    if OUTDIR == "efficiency"
        errors = collect(-0.5:0.05:0.5)
    else
        errors = collect([-0.5:0.5:0.5..., 2.0, 4.0]) # keep 2.0 to test removal
    end

    # Create the layers with the right percentages
    layers = Dict()
    for e in errors
        # Convert the error / percentage to a number of cells
        ncell = length(sp_mask)
        ntarget = (1.0 + e) * ncell
        n = round(Int, ntarget)

        # Select the threshold based on the number of cells
        ind = n < length(probsp_range) ? n : length(probsp_range)
        ind = iszero(ind) ? 1 : ind
        thresh = sort(values(probsp_range); rev=true)[ind]

        # Create the masked layer with the closest percentage
        l = convert(SDT.SDMLayer{Float64}, probsp_range .>= thresh)
        SDT.nodata!(l, 0)
        layers[e] = l
    end
    layers

    # Remove cases where layers have the same number of cells (i.e. the exact same
    # cells) This should mostly apply for species with high occupancy when the range
    # overestimation is too high. Unless we remove them, we'll end running multiple
    # times the exact same similation. It's better to deal with it later on when
    # summarizing the results.
    removed = []
    for i in length(errors):-1:2 # needs to run backwards
        # Errors to compare
        e = errors[i]
        eprev = errors[i - 1]
        # Layers
        l = layers[e]
        lprev = layers[eprev]
        # Delete if the same as previous
        if length(l) == length(lprev)
            # @info "Removing error = $e from the set of layers as it duplicates another layer"
            delete!(layers, e)
            push!(removed, e)
        end
    end
    layers

    # Add effort adjustment
    nbon_max = 500
    nbon_ref = 300
    nstep = 31
    nmaxs = [round(Int, nbon_ref * (1 + e)) for e in reverse(errors)]
    nbons = [[1, round.(Int, range(0, n; length=nstep)[Not(1)])...] for n in nmaxs]

    # Get the maximum proportion of interactions to monitor in the layer
    pmax = Dict()
    for k in keys(layers)
        idx = findall(isone, layers[k])
        nets = SIS.networks(realized)[idx]
        monitored_int = unique(reduce(vcat, interactions.(render.(Binary, nets))))
        monitored_deg = sum(in.(sp, monitored_int))
        realized_deg = degree(realized.metaweb, sp)
        # pmax[k] = monitored_deg / realized_deg
        pmax[k] = monitored_deg
    end
    pmax
    pmax_df = sort(DataFrame((id=id, offset=k, pmax=v) for (k, v) in pmax), :offset)
    return CSV.write(datadir("_xtras", "pmax", "pmax-$idp.csv"), pmax_df)
end
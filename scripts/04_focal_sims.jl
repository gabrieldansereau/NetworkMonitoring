using DrWatson
@quickactivate :NetworkMonitoring

using BiodiversityObservationNetworks:
    UncertaintySampling, BalancedAcceptance, WeightedBalancedAcceptance, SimpleRandom

# Set directory to export results
if !(@isdefined OUTDIR)
    const OUTDIR = "dev" # dev (local), focal_array or efficiency (remote)
end
mkpath(datadir(OUTDIR))

# Use job id to random results
id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
idp = lpad(id, 3, "0")

# Set number of replicates for short interactive run
if !(@isdefined NREP)
    const NREP = 3
end
# For longer, non-interactive runs, set the number of replicates with:
# julia --project -e 'const NREP = 100; include("scripts/04_focal_sims.jl")'

# Generate networks using simulations
d = DefaultParams()
seed = id * 42
nets_dict, sims_dict = generate_focal_simulation(d; seed=seed);
@unpack detected, realized, pos, metaweb, ranges = nets_dict
@unpack probranges, deg, sp, sp_range, sp_mask = sims_dict
@unpack richness_spp, richness_int, richness_pos, degree_possible, degree_realized =
    sims_dict
@unpack probsp_range, thresholds = sims_dict

# Export individual layers only for focal array simulations
SDT.SimpleSDMLayers.save(
    datadir(OUTDIR, "layer_sp_range-$idp.tiff"), [sp_range, probsp_range, sp_mask]
)
SDT.SimpleSDMLayers.save(
    datadir(OUTDIR, "layer_richness-$idp.tiff"),
    [convert(SDT.SDMLayer{Int64}, richness_spp), richness_int, richness_pos],
)
SDT.SimpleSDMLayers.save(
    datadir(OUTDIR, "layer_degree-$idp.tiff"), [degree_realized, degree_possible]
)

# Test run focal monitoring
@info "Test run"
Random.seed!(333)
monitored_sp = focal_monitoring(
    nets_dict, sp; name="test", type=[:possible], nrep=2, nbons=1:5:100
)

# Export test for convenience in plot scripts
if id == 1
    CSV.write(datadir("monitored_test.csv"), monitored_sp)
end

## Repeat focal monitoring by network types

#=

# Run for all types
@info "Network types"
Random.seed!(id * 33)
types = [:possible, :realized, :detected]
STEP = (OUTDIR == "dev" ? 20 : 5)
monitored_types = focal_monitoring(
    nets_dict, sp; name="types", type=types, nrep=NREP, nbons=1:STEP:101, combined=false
)

# Re-run for realized and detected with more sites in BON
types2 = [:realized, :detected]
STEP = (OUTDIR == "dev" ? 4000 : 1000)
monitored_types2 = focal_monitoring(
    nets_dict,
    sp;
    name="types2",
    type=types2,
    nbons=1:STEP:10_001,
    nrep=NREP,
    combined=false,
)

# Export
CSV.write(datadir(OUTDIR, "monitored_types-$idp.csv"), monitored_types)
CSV.write(datadir(OUTDIR, "monitored_types2-$idp.csv"), monitored_types2)

=#

## Repeat with 4 species with different degrees

#=

Random.seed!(id * 101)

# Get species to test
degrees = degree(metaweb.metaweb)
idx = [1, 25, 50, 70]
spp = sort(collect(degrees); by=x -> x.second, rev=true)[idx]
spp = [sp.first for sp in spp]

# Repeat focal monitoring per species
@info "Species"
STEP = (OUTDIR == "dev" ? 50 : 10)
monitored_spp = focal_monitoring(
    nets_dict,
    spp;
    name="species",
    type=[:realized],
    nrep=NREP,
    nbons=1:STEP:500,
    combined=false,
)

# Extract species ranges
speciesranges = [
    SDT.SDMLayer(occurrence(ranges)[id]; x=(0.0, d.nsites), y=(0.0, d.nsites)) for id in idx
]

# Get layer occupancy
occupancies = occupancy.(speciesranges)
monitored_spp_occ = DataFrame(; sim=id, sp=spp, rank=1:4, occ=occupancies)

# Export
CSV.write(datadir(OUTDIR, "monitored_spp-$idp.csv"), monitored_spp)
CSV.write(datadir(OUTDIR, "monitored_spp_occ-$idp.csv"), monitored_spp_occ)

=#

## Explore variations with different sampler

#=

# Run with replicates
@info "Samplers"
Random.seed!(id * 22)
samplers = [UncertaintySampling, WeightedBalancedAcceptance, SimpleRandom]
STEP = (OUTDIR == "dev" ? 50 : 5)
monitored_samplers = focal_monitoring(
    nets_dict,
    sp,
    sp_range;
    name="samplers",
    type=[:realized],
    sampler=samplers,
    nbons=1:STEP:500,
    nrep=NREP,
    combined=false,
)

# Variation with focal range mask
@info "Focal range mask"
Random.seed!(id * 23)
monitored_mask = focal_monitoring(
    nets_dict,
    sp,
    sp_mask;
    name="samplers",
    type=[:realized],
    sampler=[BalancedAcceptance, SimpleRandom],
    nbons=1:STEP:500,
    nrep=NREP,
    combined=false,
)
@rtransform!(monitored_mask, :sampler = string(:sampler))
replace!(monitored_mask.sampler, "SimpleRandom" => "SimpleRandomMask")
append!(monitored_samplers, monitored_mask; promote=true)

# Export
CSV.write(datadir(OUTDIR, "monitored_samplers-$idp.csv"), monitored_samplers)

=#

## Explore variations with optimization layers

#=

# Optimize with UncertaintySampling
@info "Optimization layers"
Random.seed!(id * 44)
optim = [richness_spp, degree_realized, probsp_range]
optimlabels = ["Species richness", "Realized interactions", "Probabilistic range"]
STEP = (OUTDIR == "dev" ? 50 : 5)
monitored_optimized = focal_monitoring(
    nets_dict,
    sp,
    optim;
    name="layers",
    type=[:realized],
    sampler=[UncertaintySampling],
    nbons=1:STEP:500,
    nrep=NREP,
    combined=false,
)
@rtransform!(monitored_optimized, :layer = optimlabels[:layer])
@select!(monitored_optimized, :set, :sp, :type, :sampler, :layer, All())

# Combine with Uncertainty Sampling on focal species layer
if !(@isdefined monitored_samplers)
    @warn "Need to run previous section for Uncertainty Sampling option"
else
    @chain begin
        monitored_samplers
        filter(:sampler => ==(UncertaintySampling), _)
        @transform(:layer = "Focal species range", :set = "layers")
        append!(monitored_optimized, _; promote=true, cols=:subset)
    end
end

# Export
CSV.write(datadir(OUTDIR, "monitored_optimized-$idp.csv"), monitored_optimized)

=#

## Range estimation

# Misestimate ranges
if OUTDIR == "efficiency"
    errors = -0.3:0.02:0.3
else
    errors = -0.2:0.05:0.2
end
threshold = thresholds[indexin([sp], probranges.species)...]
layers = Dict()
for e in errors
    l = convert(SDT.SDMLayer{Float64}, probsp_range .> threshold + e)
    SDT.nodata!(l, 0)
    layers[e] = l
end
layers

# Optimize with UncertaintySampling
@info "Range estimations"
Random.seed!(id * 832)
set = errors
optim = [layers[s] for s in set]
optimlabels = [ifelse(s < 0, "Over$s", "Under-$s") for s in set]
replace!(optimlabels, "Under-0.0" => "True-0.0")
STEP = (OUTDIR == "dev" ? 50 : 5)
monitored_estimations = focal_monitoring(
    nets_dict,
    sp,
    optim;
    name="ranges",
    type=[:realized],
    sampler=[BalancedAcceptance],
    nbons=1:STEP:500,
    nrep=NREP,
    combined=false,
)
@rtransform!(monitored_estimations, :layer = optimlabels[:layer])
@select!(monitored_estimations, :set, :sp, :type, :sampler, :layer, All())

# Export
CSV.write(datadir(OUTDIR, "monitored_estimations-$idp.csv"), monitored_estimations)
SDT.SimpleSDMLayers.save(datadir(OUTDIR, "layer_range_estimations-$idp.tiff"), optim)

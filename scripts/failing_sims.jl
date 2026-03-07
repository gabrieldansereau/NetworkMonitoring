# 2026-03-06
# Checking out what happens with missing simulations in range estimations
# It's very likely tied to occupancy and having too few presence sites for the BON

include("include.jl")

import BiodiversityObservationNetworks as BON
import NetworkMonitoring as NM

# # Preliminary setup

# Calculate occupancy
ids = collect(1:100)
occup = Dict(ids .=> zeros(length(ids)))
masked = Dict()
for i in ids
    idp = lpad(i, 2, "0")
    l = SDT.SDMLayer(datadir("efficiency", "layer_sp_range-$idp.tiff"))
    occup[i] = occupancy(l)
    masked[i] = SDT.SDMLayer(datadir("efficiency", "layer_sp_mask-$idp.tiff"))
end
occupdf = DataFrame(;
    sim=ids, occ=[occup[i] for i in ids], nocc=[length(masked[i]) for i in ids]
)
sort(occupdf, :nocc)
# all true ranges have more than 500 sites

hist(occupdf.nocc)

# All is well distributed

# do I have the info already saved in estimation results?
effs_estimations = CSV.read(datadir("efficiency_estimations.csv"), DataFrame)
sort(unique(effs_estimations, [:sim, :occ]), :occ)
# nope occupancy is the same per sim, I'll have to calculate them

## Investigate occupancies manually

# get from exported layers
errors = collect(-0.3:0.02:0.3)
estim_occup = Dict()
layers = Dict()
estimdf = DataFrame()
for i in ids
    idp = lpad(i, 2, "0")
    f = datadir("efficiency", "layer_range_estimations-$idp.tiff")
    isfile(f) || continue
    layers[i] = Dict()
    for (j, e) in enumerate(errors)
        l = SDT.SDMLayer(f; bandnumber=j)
        o = length(l)
        estim_occup["$idp-$j"] = o
        layers[i][e] = l
        push!(estimdf, (; sim=i, n=j, err=e, occ=o, layer=l))
    end
end
estim_occup
layers
estimdf

probs = @rsubset(estimdf, :occ < 500)
@chain estimdf begin
    @rsubset(:occ < 500)
    groupby(:sim)
    combine(nrow => :ncases, :occ => minimum => :min)
    sort(:ncases; rev=true)
    # sort(:min)
end

# Failed layers are not included here as they don't export their files
# Layer 39 will be a good test case with a minimum of 496 sites
# Layer 59 is the complete opposite with a minimum of 1... --> timeout in the first array, ran after
# Layer 83 timed out too, but has a few more sites
# Layer 100 is the successful one with the most cases < 500
# Maybe other layers fail because they have minimums of 0?

ltest0 = layers[39][0.0]
ltest1 = layers[39][0.3]
ltest2 = layers[83][0.3]
ltest3 = layers[59][0.3]

lnull = similar(ltest3)
SDT.nodata!(lnull, 0.0)

## Explore BONs

Random.seed!(42);
samp = BON.BalancedAcceptance
@time b0 = BON.sample(samp(500), ltest0) # test case
@time b1 = BON.sample(samp(500), ltest1) # no issues with 496
@time b2 = BON.sample(samp(500), ltest2) # 33 is fine
@time b3 = BON.sample(samp(500), ltest3) # 1 takes 7 sec but works
# DO NOT RUN! @time bn = BON.sample(BON.BalancedAcceptance(500), lnull) # takes forever
# so the issue is when there are no sites in the layer, which is simple to verify

b0.nodes |> unique # all unique
b1.nodes |> unique # not sampling all possible sites
b2.nodes |> unique # so restricted it actually does sample all possible sites
b3.nodes |> unique # same as 33
unique(b3) # same
b3.nodes # kinda silly, should at least call unique to shorten

NM.coordinates(b0) |> unique
NM.coordinates(b1) |> unique
NM.coordinates(b2) |> unique
NM.coordinates(b3) |> unique
# same as using nodes

## Repeat with different samplers

# UncertaintySampling: super fast with all (not actuallyu used with masked layers)
Random.seed!(42);
samp = BON.UncertaintySampling
@time b0 = BON.sample(samp(500), ltest0)
@time b1 = BON.sample(samp(500), ltest1)
@time b2 = BON.sample(samp(500), ltest2)
@time b3 = BON.sample(samp(500), ltest3) # actually faster

# WeightedBalancedAcceptance: DO NOT TRY, not converging either?
Random.seed!(42);
samp = BON.WeightedBalancedAcceptance
@time b0 = BON.sample(samp(500), ltest0) # not converging either?
@time b1 = BON.sample(samp(500), ltest1)
@time b2 = BON.sample(samp(500), ltest2)
@time b3 = BON.sample(samp(500), ltest3)

# SimpleRandom
Random.seed!(42);
samp = BON.SimpleRandom
@time b0 = BON.sample(samp(500), ltest0) # test case
@time b1 = BON.sample(samp(500), ltest1) # no issues with 496
@time b2 = BON.sample(samp(500), ltest2) # 33 is fine
@time b3 = BON.sample(samp(500), ltest3)

## Check failed sims

# Alright so let's see about those failed sims

fails = setdiff(ids, unique(estimdf.sim)) # [24, 44, 74]

errors = -0.3:0.02:0.3
failed_thresholds = Dict()
failed_layers = Dict()
for id in fails
    # Generate networks using simulations
    d = NM.DefaultParams()
    seed = id * 42
    nets_dict, sims_dict = NM.generate_focal_simulation(d; seed=seed)
    @unpack sp, probranges, thresholds, probsp_range = sims_dict

    # Misestimate ranges
    t = thresholds[indexin([sp], probranges.species)...]
    fl = [convert(SDT.SDMLayer{Float64}, probsp_range .> t + e) for e in errors]
    SDT.nodata!.(fl, 0)

    # Extract
    failed_thresholds[id] = t
    failed_layers[id] = fl
end

collect(values(failed_thresholds)) # weird how they are the same...

collect(values(failed_layers))

failed_ends = Dict(k => last(length.(v), 5) for (k, v) in failed_layers)

# 44 and 74 have zeros, which explains why they fail
# 24 does not, but it did succeed once in 190 min with the original seed
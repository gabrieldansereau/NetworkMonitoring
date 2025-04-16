using NetworkMonitoring

# Random.seed!(42)

## Generate network & abundances

# Define all required variables
# set_params = true
if !(@isdefined set_params) || set_params
    @info "Setting parameters to default values"
    ns = 100
    nsites = 100
    C_exp = 0.2
    ra_sigma = 1.2
    ra_scaling = 50.0
    energy_NFL = 10_000
    H_nlm = 0.5
    nbon = 50
end;

# Generate metaweb using Niche Model
# Random.seed!(42);
metaweb = generate(SIS.NicheModel(ns, C_exp))

# Generate autocorrelated ranges
ranges = generate(AutocorrelatedRange(dims=(nsites, nsites)),ns)

# Generate realized abundances
ra = generate(NormalizedLogNormal(σ=ra_sigma), metaweb)

# Generate possible local networks
pos = possible(metaweb, ranges)

# Generate detectable metaweb given abundances
detectable = detectability(RelativeAbundanceScaled(ra_scaling), metaweb, ra)

# Generate possible spatial metaweb (given cooccurrence)
realized = possible(metaweb, ranges)
# Make network realizable given Neutrally forbidden links
realizable!(NeutrallyForbiddenLinks(energy_NFL), realized, ra)
# Realize it
realize!(realized)

# Detect it
detected = deepcopy(realized)
detect!(detected, detectable)

## Extract network info - Manually from scale objects

# Build detected metaweb
detected_metaweb = metawebify(detected)

# Build realized metaweb
realized_metaweb = metawebify(realized)

# Build possible metaweb
possible_metaweb = metawebify(pos)

# Build binary metaweb
binary_metaweb = Matrix{Int}(adjacency(metaweb))

## Extract network measures

# Extract detected link matrix
detected_links = extract(links, detected)

# Extract realized link matrix
realized_links = extract(links, realized)

# Extract possible link matrix
possible_links = extract(links, pos)

# Extract richess
ranges_richness = sum(occurrence(ranges))

## Add BON

# Generate uncertainty layer
uncertainty = SDT.SDMLayer(
    MidpointDisplacement(H_nlm), (nsites, nsites);
    x=(0.0, nsites), y=(0.0, nsites)
)

# Use BalancedAcceptance sampling
bon = BON.sample(BON.BalancedAcceptance(nbon), uncertainty)
boncoords = coordinates(bon)
bonidx = [get_grid_coordinate_by_latlon(uncertainty, bc...) for bc in boncoords]
xs = getfield.(bonidx, 1)
ys = getfield.(bonidx, 2)

## Extract infos from BON

# Extract sites
sites = DataFrame(x=xs, y=ys)
sites.coords = CartesianIndex.(xs, ys)

# Extract richness/species info
@transform!(sites, :richness = ranges_richness[sites.coords])
# Extract species observed from ranges
ranges_mat = Array{Union{Nothing, Symbol}}(nothing, ns, ns, nsites)
for i in 1:ns
    ranges_mat[findall(!iszero, ranges.occurrence[i].range_map), i] .= ranges.species[i]
end
ranges_mat
# Extract species observed at sites
species_mat = [unique(filter(!isnothing, ranges_mat[x, y, :])) for x in 1:ns, y in 1:ns]
@transform!(sites, :species_set = [species_mat[r.x, r.y] for r in eachrow(sites)])
# List observed species
observed_species = Vector{Symbol}(unique(reduce(vcat, sites.species_set)))
prop_observed_sp = length(observed_species) / ns
# @info "Proportion of species observed at sampled sites:
#     $(round(prop_observed_sp; sigdigits=3))"

# Extract links/interaction info
@transform!(
    sites,
    :links_detected = detected_links[sites.coords],
    :links_realized = realized_links[sites.coords],
    :links_possible = possible_links[sites.coords],
)
# Extract local networks at sites
@transform!(
    sites,
    :detected = render.(Binary, detected.scale.network[sites.coords]),
    :realized = render.(Binary, realized.scale.network[sites.coords]),
    :possible = render.(Binary, pos.scale.network[sites.coords]),
)
# Extract interactions at sites
@transform!(
    sites,
    :int_detected = interactions.(sites.detected),
    :int_realized = interactions.(sites.realized),
    :int_possible = interactions.(sites.possible),
)

# List unique interactions
sampled_int_detected = unique(reduce(vcat, sites.int_detected))
sampled_int_realized = unique(reduce(vcat, sites.int_realized))
sampled_int_possible = unique(reduce(vcat, sites.int_possible))

# List all interactions
all_int_detected = interactions(render(Binary, detected.metaweb))
all_int_realized = interactions(render(Binary, realized.metaweb))
all_int_possible = interactions(render(Binary, pos.metaweb))

# List proportions
prop_detected_int = length(sampled_int_detected) / length(all_int_detected)
prop_realized_int = length(sampled_int_realized) / length(all_int_realized)
prop_possible_int = length(sampled_int_possible) / length(all_int_possible)

# Info message
#=
@info "Proportion of sampled interactions based on detected interactions:
    $(round(prop_detected_int; sigdigits=3))"
@info "Proportion of sampled interactions based on realized interactions:
    $(round(prop_realized_int; sigdigits=3))"
@info "Proportion of sampled interactions based on possible interactions:
    $(round(prop_possible_int; sigdigits=3))"
 =#
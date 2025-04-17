using NetworkMonitoring

# Random.seed!(42)

## Generate network & abundances

# Define all required variables
# set_params = true
if !(@isdefined set_params) || set_params
    @info "Setting parameters to default values"
    ns = 75
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
possible_richness = extract(sum, ranges)

## Add BON

# Generate uncertainty layer
uncertainty = SDT.SDMLayer(
    MidpointDisplacement(H_nlm), (nsites, nsites);
    x=(0.0, nsites), y=(0.0, nsites)
)

# Use BalancedAcceptance sampling
bon = BON.sample(BON.BalancedAcceptance(nbon), uncertainty)

## Extract infos from BON

# List observed species
sp_monitored = monitor(x -> findall(isone, x), ranges, bon; makeunique=true)
sp_monitored = NetworkMonitoring._getspecies(sp_monitored, ranges)

# List monitored interactions
int_detected = monitor(x -> interactions(render(Binary, x)), detected, bon; makeunique=true)
int_realized = monitor(x -> interactions(render(Binary, x)), realized, bon; makeunique=true)
int_possible = monitor(x -> interactions(render(Binary, x)), pos, bon; makeunique=true)

# List all interactions
all_detected = interactions(render(Binary, detected.metaweb))
all_realized = interactions(render(Binary, realized.metaweb))
all_possible = interactions(render(Binary, pos.metaweb))

# List proportions to compare
prop_detected_int = length(int_detected) / sum(detected_metaweb)
prop_realized_int = length(int_realized) / sum(realized_metaweb)
prop_possible_int = length(int_possible) / sum(possible_metaweb)
prop_monitored_sp = length(sp_monitored) / ns

# Info message
#=
@info "Proportion of monitored interactions based on detected interactions:
    $(round(prop_detected_int; sigdigits=3))"
@info "Proportion of monitored interactions based on realized interactions:
    $(round(prop_realized_int; sigdigits=3))"
@info "Proportion of monitored interactions based on possible interactions:
    $(round(prop_possible_int; sigdigits=3))"
@info "Proportion of monitored species over all species:
    $(round(prop_monitored_sp; sigdigits=3))"
 =#
using CairoMakie
using DataFramesMeta
using GeoMakie
using NeutralLandscapes
using Random
using SparseArrays
using SpeciesInteractionSamplers
using Statistics
import BiodiversityObservationNetworks as BON
import SpeciesDistributionToolkit as SDT
import SpeciesInteractionNetworks as SIN

Random.seed!(42)

## Generate network & abundances

# Generate metaweb using Niche Model
ns = 100
Random.seed!(42); metaweb = generate(NicheModel(ns,0.2))

# Generate autocorrelated ranges
nsites = 100
ranges = generate(AutocorrelatedRange(dims=(nsites, nsites)),ns)
ranges.occurrence[2].range_map

# Generate realized abundances
ra = generate(NormalizedLogNormal(σ=1.2), metaweb)
ra.abundance

# Generate possible local networks
pos = possible(metaweb, ranges)
pos.metaweb

# Generate detectable metaweb given abundances
detectable = detectability(RelativeAbundanceScaled(), metaweb, ra)

# Generate possible spatial metaweb (given cooccurrence)
realized = possible(metaweb, ranges)
# Make network realizable given Neutrally forbidden links
realizable!(NeutrallyForbiddenLinks(10000), realized, ra)
# Realize it
realize!(realized)

# Detect it
detected = deepcopy(realized)
detect!(detected, detectable)

## Extract network info - SIS.jl's metaweb objects

# Investigate output
Matrix{Int}(detected.metaweb.edges.edges) == Matrix{Int}(realized.metaweb.edges.edges)
adjacency(detected.metaweb) == adjacency(realized.metaweb)

# Extract metawebs
detected_metaweb1 = Matrix{Int}(adjacency(detected.metaweb) .>= 1)
realized_metaweb1 = Matrix{Int}(adjacency(realized.metaweb) .>= 1)
possible_metaweb1 = Matrix{Int}(adjacency(pos.metaweb) .>= 1)
binary_metaweb1 = Matrix{Int}(adjacency(metaweb))

# Check differences
detected_metaweb1 == binary_metaweb1 # different
detected_metaweb1 == possible_metaweb1 # different
detected_metaweb1 == realized_metaweb1 # same? 🤨
sum(detected_metaweb1)
sum(realized_metaweb1)
sum(possible_metaweb1)
sum(binary_metaweb1)
heatmap(binary_metaweb1 - detected_metaweb1)

# Is this worth anything?
sum(detected_metaweb1) == SIN.links(detected.metaweb) # yep
sum(realized_metaweb1) == SIN.links(realized.metaweb) # yep
sum(possible_metaweb1) == SIN.links(pos.metaweb) # yep
# So just extra steps not to gain much essentially
# Except to notice difference between detected & realized?

## Extract network info - Manually from scale objects

# Create function to extract metaweb from scale objects
function metawebify(m::T; binary=true) where T <: Metaweb
    m_adj = convert.(Matrix{Int}, adjacency(m.scale))
    m_acc = accumulate(+, vec(m_adj))
    m_end = m_acc[end]
    if binary
        m_end = Matrix{Int}(m_end .>= 1)
    end
    return m_end
end

# Build detected metaweb
detected_metaweb = metawebify(detected)
detected_metaweb == detected_metaweb1 # false? 🤨
(sum(detected_metaweb), sum(detected_metaweb1)) # why so few?
# same issue as identified previously:
#   somehow detected_metaweb1 == realized_metaweb, but should it?

# Build realized metaweb
realized_metaweb = metawebify(realized)
realized_metaweb == realized_metaweb1 # true 🥳
(sum(realized_metaweb), sum(realized_metaweb1))

# Build possible metaweb
possible_metaweb = metawebify(pos)
possible_metaweb == possible_metaweb1 # true 🥳
(sum(possible_metaweb), sum(possible_metaweb1))

# Build binary metaweb
binary_metaweb = binary_metaweb1
sum(binary_metaweb)

## Extract network measures

# Extract spatial link matrix
detected_links = zeros(Int, size(detected))
for i in 1:size(detected)[1], j in 1:size(detected)[2]
    detected_links[i,j] = SIN.links(detected.scale.network[i,j])
end
detected_links
heatmap(detected_links)

# Extract spatial link matrix
realized_links = zeros(Int, size(realized))
for i in 1:size(realized)[1], j in 1:size(realized)[2]
    realized_links[i,j] = SIN.links(realized.scale.network[i,j])
end
realized_links
heatmap(realized_links; axis=(; aspect=1))
heatmap(realized_links - detected_links; axis=(; aspect=1))

## Add BON

# Generate uncertainty layer
uncertainty = begin Random.seed!(42); rand(MidpointDisplacement(), (100, 100)) end
heatmap(uncertainty; axis=(; aspect=1))

# Generate BON
# bon = BON.sample(BON.AdaptiveHotspot(), uncertainty)
# xs = [n.coordinate[1] for n in bon.nodes]
# ys = [n.coordinate[2] for n in bon.nodes]
# begin
#     f, ax, p = heatmap(uncertainty; axis=(; aspect=1))
#     scatter!(100xs, 100ys; color="red")
#     Colorbar(f[1,end+1], p)
#     f
# end # ???

# Use Simple Random sampling instead
bon = BON.sample(BON.SimpleRandom(50), uncertainty)
xs = round.(Int, 100*[n.coordinate[1] for n in bon.nodes])
ys = round.(Int, 100*[n.coordinate[2] for n in bon.nodes])
begin
    f, ax, p = heatmap(uncertainty; axis=(; aspect=1))
    scatter!(xs, ys; color="red")
    Colorbar(f[1,end+1], p)
    f
end

# View sampling over richness
ranges_richness = sum(occurrence(ranges))
begin
    f, ax, p = heatmap(ranges_richness; axis=(; aspect=1), colormap=:cividis)
    scatter!(xs, ys; color="red")
    Colorbar(f[1,end+1], p)
    f
end

# View sampling over interactions
begin
    f, ax, p = heatmap(realized_links; axis=(; aspect=1), colormap=:cividis)
    scatter!(xs, ys; color="red")
    Colorbar(f[1,end+1], p)
    f
end

## Extract infos from BON

# Extract sites
sites = DataFrame(x=xs, y=ys)
filter!(:y => !iszero, sites) # remove site with zero as coordinate (for now)

# Extract richness/species info
@transform!(sites, :richness = [ranges_richness[x, y] for (x,y) in zip(sites.x, sites.y)])
# Extract species observed from ranges
ranges_mat = zeros(Int, ns, ns, nsites)
for i in 1:ns
    ranges_mat[:, :, i] = deepcopy(ranges.occurrence[i].range_map)
    ranges_mat[:, :, i] = replace(ranges_mat[:, :, i], 1 => i)
end
ranges_mat
# Extract species observed at sites
species_mat = [unique(filter(!iszero, ranges_mat[x, y, :])) for x in 1:ns, y in 1:ns]
@transform!(sites, :species_set = [species_mat[r.x, r.y] for r in eachrow(sites)])
# List observed species
observed_species = unique(reduce(vcat, sites.species_set))
proportion_observed_species = length(observed_species) / ns
@info "Proportion of species observed at sampled sites:
    $(round(proportion_observed_species; sigdigits=3))"

# Extract links/interaction info
@transform!(sites, :links_realized = [realized_links[r.x, r.y] for r in eachrow(sites)])
@transform!(sites, :links_detected = [detected_links[r.x, r.y] for r in eachrow(sites)])
# Extract local networks at sites
@transform!(
    sites,
    :realized = [SIN.render(SIN.Binary, realized.scale.network[r.x, r.y]) for r in eachrow(sites)],
    :detected = [SIN.render(SIN.Binary, detected.scale.network[r.x, r.y]) for r in eachrow(sites)],
)
# Extract interactions at sites
@transform!(sites, :int_realized = [SIN.interactions(r.realized) for r in eachrow(sites)])
@transform!(sites, :int_detected = [SIN.interactions(r.detected) for r in eachrow(sites)])

# List unique interactions
observed_int = unique(reduce(vcat, sites.int_realized))
detected_int = unique(reduce(vcat, sites.int_detected))
realized_int = SIN.interactions(SIN.render(SIN.Binary, realized.metaweb))

# List proportions
proportion_observed_int = length(observed_int) / length(realized_int)
@info "Proportion of interactions observed at sampled sites:
    $(round(proportion_observed_int; sigdigits=3))"
proportion_detected_int = length(detected_int) / length(realized_int)
@info "Proportion of interactions detected at sampled sites:
    $(round(proportion_detected_int; sigdigits=3))"


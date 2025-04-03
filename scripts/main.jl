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

# Investigate output
Matrix{Int}(detected.metaweb.edges.edges) == Matrix{Int}(realized.metaweb.edges.edges)
adjacency(detected.metaweb) == adjacency(realized.metaweb)
detected_metaweb = Matrix{Int}(adjacency(detected.metaweb) .>= 1)
realized_metaweb = Matrix{Int}(adjacency(realized.metaweb) .>= 1)
binary_metaweb = Matrix{Int}(adjacency(metaweb))
detected_metaweb == binary_metaweb # different
detected_metaweb == realized_metaweb # same 🤨
sum(detected_metaweb)
sum(binary_metaweb)
heatmap(binary_metaweb - detected_metaweb)

# Build detected metaweb
d_all = convert.(Matrix{Int}, adjacency(detected.scale))
d_acc = accumulate(+, vec(d_all))
detected_metaweb = Matrix{Int}(d_acc[end] .>= 1)
sum(detected_metaweb)
sum(binary_metaweb)
heatmap(binary_metaweb - detected_metaweb; axis=(; aspect=1))

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

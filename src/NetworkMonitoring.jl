module NetworkMonitoring

using Reexport

@reexport begin
    using AlgebraOfGraphics
    using CairoMakie
    using CSV
    using DrWatson
    using DataFramesMeta
    using GeoMakie
    using NeutralLandscapes
    using ProgressMeter
    using Random
    using SparseArrays
    using SpeciesInteractionNetworks
    using SpeciesInteractionSamplers
    using Statistics

    import BiodiversityObservationNetworks as BON
    import SpeciesDistributionToolkit as SDT
    import SpeciesDistributionToolkit.SimpleSDMLayers as SSL
    import SpeciesInteractionSamplers as SIS

    using BiodiversityObservationNetworks: GI.coordinates
    using SpeciesDistributionToolkit: SimpleSDMLayers.__get_grid_coordinate_by_latlon as get_grid_coordinate_by_latlon
    using SpeciesInteractionNetworks: SpeciesInteractionNetworks as SIN, interactions
end

include("metaweb.jl")
include("monitor.jl")
include("plots.jl")
include("ranges.jl")
include("utils.jl")

# Trick LanguageServer into cooperating in script files
@static if false
    include("../scripts/main.jl")
    include("../scripts/01_param_search.jl")
    include("../scripts/02_figures.jl")
end

export BON, SDT, SSL, SIS, SIN
export metawebify, extract, heatmapcb, monitor

end # module NetworkMonitoring

module NetworkMonitoring

using Reexport

@reexport begin
    # Load packages
    using AlgebraOfGraphics
    using CairoMakie
    using CSV
    using DataFramesMeta
    using DrWatson
    using Distributions
    using GeoMakie
    using NeutralLandscapes
    using ProgressMeter
    using Random
    using SparseArrays
    using SpeciesInteractionNetworks
    using SpeciesInteractionSamplers
    using Statistics

    # Import packages with shorter names for convenience
    import BiodiversityObservationNetworks as BON
    import SpeciesDistributionToolkit as SDT
    import SpeciesDistributionToolkit.SimpleSDMLayers as SSL
    import SpeciesInteractionSamplers as SIS

    # Load specific functions
    using BiodiversityObservationNetworks: GI.coordinates
    using SpeciesDistributionToolkit:
        SimpleSDMLayers.__get_grid_coordinate_by_latlon as get_grid_coordinate_by_latlon
    using SpeciesInteractionNetworks: SpeciesInteractionNetworks as SIN, interactions

    # Import functions to extend
    import SpeciesInteractionSamplers: generate
end

include("metaweb.jl")
include("monitor.jl")
include("plots.jl")
include("ranges.jl")
include("runsim.jl")
include("utils.jl")

# Trick LanguageServer into cooperating in script files
@static if false
    # include("../scripts/file.jl")
end

export BON, SDT, SSL, SIS, SIN
export metawebify, extract, heatmapcb, monitor
export DefaultParams, generate_networks, generate_bon, evaluate_monitoring, runsim
export AutocorrelatedProbabilisticRange, generate

end # module NetworkMonitoring

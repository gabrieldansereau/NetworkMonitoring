"""
Extract biodiversity information from monitored sites using a BON.
"""
function monitor(m::T, bon::BON.BiodiversityObservationNetwork) where T <: Metaweb
    # Extract sites
    coords = coordinates(bon)

    # Extract site indices
    _x, _y = size(m)
    layer = SDT.SDMLayer(
        zeros(Bool, (_x, _y));
        x=(0.0, _x), y=(0.0, _y)
    )
    idx = [CartesianIndex(get_grid_coordinate_by_latlon(layer, c...)) for c in coords]

    # Extract network info
    nets = SIS.networks(m)[idx]

    return nets
end

function monitor(
    f::Function, m::T, bon::BON.BiodiversityObservationNetwork; makeunique::Bool=false
) where T <: Metaweb
    monitored = f.(monitor(m, bon))
    if makeunique
        return unique(reduce(vcat, monitored))
    else
        return monitored
    end
end

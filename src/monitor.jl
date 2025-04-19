"""
    monitor(m::T, bon::BON.BiodiversityObservationNetwork)  where T <: Metaweb
    monitor(r::T, bon::BON.BiodiversityObservationNetwork)  where T <: Occurrence

Extract biodiversity information from monitored sites using a BON.
"""
function monitor(m::T, bon::BON.BiodiversityObservationNetwork) where T <: Metaweb
    # Extract sites
    coords = coordinates(bon)

    # Extract site indices
    _x, _y = size(m.scale)
    layer = SDT.SDMLayer(
        zeros(Bool, (_x, _y));
        x=(0.0, _x), y=(0.0, _y)
    )
    idx = [CartesianIndex(get_grid_coordinate_by_latlon(layer, c...)) for c in coords]

    # Extract network info
    nets = SIS.networks(m)[idx]

    return nets
end
function monitor(r::T, bon::BON.BiodiversityObservationNetwork) where T <: Occurrence
    # Extract sites
    coords = coordinates(bon)

    # Extract site indices
    _x, _y = size(r)
    layer = SDT.SDMLayer(
        zeros(Bool, (_x, _y));
        x=(0.0, _x), y=(0.0, _y)
    )
    idx = [CartesianIndex(get_grid_coordinate_by_latlon(layer, c...)) for c in coords]

    # Extract occurrence info
    s = stack(occurrence(r))
    monitored = [s[id, :] for id in idx]

    return monitored
end


"""
    monitor(
        f::Function, r::T, bon::BON.BiodiversityObservationNetwork; makeunique::Bool=false
    ) where T <: Union{Occurrence, Metaweb}

Extract biodiversity information from monitored sites using a BON and directly
apply function `f` over extracted elements. The default value of
`makeunique=false` keeps separate information for monitored sites. Set
`makeunique=true` to reduce to unique values, for instance to extract all
monitored species or interactions across monitored sites.
"""
function monitor(
    f::Function, x::T, bon::BON.BiodiversityObservationNetwork; makeunique::Bool=false
) where T <: Union{Occurrence, Metaweb}
    monitored = f.(monitor(x, bon))
    if makeunique
        return unique(reduce(vcat, monitored))
    else
        return monitored
    end
end

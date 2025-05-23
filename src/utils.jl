function _getspecies(vals, sp_list::Vector)
    spdict = Dict(i => sp for (i, sp) in enumerate(sp_list))
    spp = getindex.(Ref(spdict), vals)
    return spp
end
_getspecies(vals, r::T) where {T<:Occurrence} = _getspecies(vals, r.species)
_getspecies(vals, m::T) where {T<:Metaweb} = _getspecies(vals, m.species.names)
"""
Extract information *site-wise* from an `Occurrence{Range}` object. Note that
some functions may perform efficiently on the output of `occurrence`, e.g.,
`sum(occurrence(ranges))`.

To extract information *species-wise* instead, use broadcasting with
`occurrence`, e.g. `sum.(occurrence(ranges))`.
"""
function extract(f::Function, r::T) where T <: Occurrence
    s = stack(occurrence(r))
    return dropdims(mapslices(f, s; dims=3), dims=3)
end

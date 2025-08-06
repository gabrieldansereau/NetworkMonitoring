"""
Extract information *site-wise* from an `Occurrence{Range}` object. Note that
some functions may perform efficiently on the output of `occurrence`, e.g.,
`sum(occurrence(ranges))`.

To extract information *species-wise* instead, use broadcasting with
`occurrence`, e.g. `sum.(occurrence(ranges))`.
"""
function extract(f::Function, r::T) where {T<:Occurrence}
    s = stack(occurrence(r))
    return dropdims(mapslices(f, s; dims=3); dims=3)
end

# ========================================
#
# Generators
#
# ========================================
"""
    AutocorrelatedProbabilisticRange <: RangeGenerator

A [`RangeGenerator`](@ref) that uses NeutralLandscapes.jl's `DiamondSquare`
landscape generator to create autocorrelated rasters with autocorrelatation
parameter ranging between 0 and 1, where increasing values mean increasing
autocorrelated. The autocorrelated raster is _not_ thresholded by default (as in
an `AutocorrelatedRange`), but the `threshold` field is kept for compatibility
when using random seeds. Setting `binary = true` will use the threshold values
and is equivalent to using `AutocorrelatedRange`.
"""
@kwdef struct AutocorrelatedProbabilisticRange <: RangeGenerator
    autocorrelation = 0.85
    dims = (50, 50)
    threshold = Distributions.Beta(10, 10)
    binary = false
end

"""
    generate(ar::AutocorrelatedProbabilisticRange)

Generates a [`Range`](@ref) using the [`AutocorrelatedProbabilisticRange`](@ref)
[`RangeGenerator`](@ref).
"""
function generate(ar::AutocorrelatedProbabilisticRange)
    H, sz, thres, bin = ar.autocorrelation, ar.dims, rand(ar.threshold), ar.binary

    range_mat = rand(DiamondSquare(H), sz)

    if bin
        range_mat[findall(x -> x < thres, range_mat)] .= 0
        range_mat[findall(!iszero, range_mat)] .= 1
    end

    return Range(range_mat)
end

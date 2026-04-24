## Related functions

# Layer occupancy
occupancy(l; f=isone) = length(findall(f, l)) / length(l.grid)

# NDI
ndi(x, y) = (x - y) / (x + y)

## Efficiency

# Saturation equation
saturation(a, pmax=1.0) = (x) -> pmax .* x ./ (a .+ x)

# Grid search for best fit
function efficiency_gridsearch(
    x, y, pmax=1.0; A=LinRange(-5.0, 15.0, 10_000), rmse=false, f=identity
)
    A = f.(A)
    err = zeros(length(A))
    for i in eachindex(A)
        f = saturation(A[i], pmax)
        err[i] = sqrt(mean((y .- f(x)) .^ 2.0))
    end
    if rmse
        _err, _ind = findmin(err)
        # _rmse = _err / length(x)
        _rmse = _err
        return (; efficiency=A[_ind], rmse=_rmse)
    else
        return A[last(findmin(err))]
    end
end

# Integral of efficiency curve
efficiency_integral(a, k=10_000, pmax=1.0) = (pmax * (a * log(a) - a * log(a + k) + k))

# Difference between two curves
function efficiency_difference(a1, a2; k=10_000, pmax=1.0)
    return efficiency_integral(a1, k, pmax) - efficiency_integral(a2, k, pmax)
end

# Value n at which the curve reaches proportion p
efficiency_n_at_p(a, p, pmax=1.0) = p * a / (pmax - p)

# Complete efficiency worklow
function efficiency(
    x,
    y;
    pmax=1.0,
    option=:integral,
    k=10_000,
    p::Float64=0.95,
    n::Int=maximum(x),
    rmse=false,
    kw...,
)
    opts = [
        :integral,
        :integral_at_n,
        :n_at_p,
        :n_at_pmax0,
        :n_at_pmax1,
        :n_at_pmax2,
        :n_at_pmax3,
        :n_at_pmax4,
        :n_at_pmax5,
        :p_at_n,
        :a,
    ]
    @assert option in opts || error("possible values for keyword option are $opts")
    if option == :n_at_pmax0
        pmax = 1.0
        option = :n_at_pmax1
    end
    a, _rmse = efficiency_gridsearch(x, y, pmax; rmse=true, kw...) # rmse=true on purpose
    if option == :integral
        ei = efficiency_integral(a, k, pmax)
    elseif option == :integral_at_n
        ei = efficiency_integral(a, n, pmax)
    elseif option == :n_at_p
        if p > pmax
            @warn "Requested p > pmax (pmax=$pmax). Defaulting to p=0.99*pmax"
            p = 0.99pmax
        end
        ei = efficiency_n_at_p(a, p, pmax)
    elseif option == :n_at_pmax1
        ei = efficiency_n_at_p(a, p * pmax, pmax)
    elseif option == :n_at_pmax2
        ei = efficiency_n_at_p(a, p * pmax, pmax) / pmax
    elseif option == :n_at_pmax3
        ei = (1 + (1 - pmax)) * efficiency_n_at_p(a, p * pmax, pmax)
    elseif option == :n_at_pmax4
        # reverse-engineer of what we would get with pmax=false
        a_under_pmax = a
        a_under = efficiency_gridsearch(x, y, 1.0; kw...) # rmse=true on purpose
        ei = a_under / a_under_pmax * efficiency_n_at_p(a_under_pmax, p * pmax, pmax)
        # which is the same as option :n_at_p (when pmax=false to get a = a_under)
        # and where we compute ei = efficiency_n_at_p(a, p, 1.0)
    elseif option == :n_at_pmax5
        # same as :n_at_pmax4 but only correcting when a_under_pmax < a_under (more efficient)
        a_under_pmax = a
        a_under = efficiency_gridsearch(x, y, 1.0; kw...) # rmse=true on purpose
        if a_under_pmax < a_under
            ei = a_under / a_under_pmax * efficiency_n_at_p(a_under_pmax, p * pmax, pmax)
        else
            ei = efficiency_n_at_p(a_under_pmax, p * pmax, pmax)
        end
        # which is the same as option :n_at_p (when pmax=false to get a = a_under)
        # and where we compute ei = efficiency_n_at_p(a, p, 1.0)
    elseif option == :p_at_n
        ei = saturation(a, pmax)(n)
    elseif option == :a
        ei = a
    end
    if rmse
        return (; ei=ei, rmse=_rmse)
    else
        return ei
    end
end

# Compare efficiencies within simulations
function comparewithin(
    effs_combined, set; to=nothing, labels=Dict(set .=> set), f=(x, y) -> -(x, y), vref=1.0
)
    # Make sure all labels are defined
    all_labels = Dict(s in keys(labels) ? s => labels[s] : s => s for s in set)
    # Select variables in set & prepare for comparison
    gpvars = [:sim, :set, :occ]
    effs_selected = @chain effs_combined begin
        # Select main columns
        select(gpvars, :variable, :eff, :eff_low, :eff_upp)
        # Select variables
        @rsubset :variable in set
        # Rename
        @rtransform :variable = all_labels[:variable]
        # Nest efficiencies in single column
        @rtransform :eff = (eff=:eff, low=:eff_low, upp=:eff_upp)
        select(gpvars, :variable, :eff)
    end
    # Define combinations to compare
    if isnothing(to)
        # By default we select all unique combinations
        comps = collect(combinations(set, 2))
    else
        # With the `to` keyword we can select a specific variable for all comparisons
        comps = []
        if !(to isa Vector)
            to = [to]
        end
        for t in to
            @assert t in set "variables selected with `to` must be in `set`"
            allcomps = vec([[v2, v1] for v1 in set, v2 in set])
            comp = filter(x -> first(x) == t, allcomps)
            filter!(x -> x[1] != x[2], comp)
            push!(comps, comp)
        end
        comps = reduce(vcat, comps)
    end
    # Compare efficiencies
    comps_df = DataFrame()
    for (c1, c2) in comps
        # Select variables for comparison
        l1 = all_labels[c1]
        l2 = all_labels[c2]
        comp_view = @rsubset(effs_selected, :variable in [l1, l2]; view=true)
        # Unstack to compare column-wise
        comp = unstack(comp_view, :variable, :eff)
        rename!(comp, Dict(l1 => :v1, l2 => :v2))
        # Create empty columns to hold results
        comp.variable .= "Δ$(l1)_$(l2)"
        comp.value = Vector{Union{Missing,Float64}}(missing, nrow(comp))
        comp.overlap = Vector{Union{Missing,Bool}}(missing, nrow(comp))
        comp.overlap_sign = Vector{Union{Missing,String}}(missing, nrow(comp))
        # Loop to compare individual results
        for r in eachrow(comp)
            eff1 = r.v1
            eff2 = r.v2
            if !ismissing(eff1) && !ismissing(eff2)
                # Actual comparison
                r.value = f(eff1.eff, eff2.eff)
                # Overlap of confidence intervals
                r.overlap = eff1.low <= eff2.upp && eff2.low <= eff1.upp
                r.overlap_sign =
                    r.overlap ? "overlap" : (r.value <= vref ? "negative" : "positive")
            end
        end
        # Simplify & export
        select!(comp, :sim, :set, :occ, :variable, :value, :overlap_sign => :overlap)
        append!(comps_df, comp)
    end
    # Remove missing comparisons
    dropmissing!(comps_df)
    disallowmissing!(comps_df)
    return comps_df
end

# Flip some comparison values
function flipthatcomp!(df, toflip; f=(x) -> -x)
    for comp in toflip
        # Arrange new comparison
        v1, v2 = split(replace(comp, "Δ" => ""), "_")
        newcomp = "Δ$(v2)_$(v1)"
        # Get current indices
        inds = findall(==(comp), df.variable)
        # Update
        new = @view df[inds, :]
        @rtransform!(
            new,
            :variable = newcomp,
            :value = f(:value),
            :overlap = replace(:overlap, "positive" => "negative", "negative" => "positive")
        )
    end
    return df
end
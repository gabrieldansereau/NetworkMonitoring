# Saturation equation
saturation(a) = (x) -> x ./ (exp(a) .+ x)

# Efficiency grid search
function efficiency_gridsearch(x, y; A=LinRange(-5.0, 15.0, 10_000), rmse=false)
    err = zeros(length(A))
    for i in eachindex(A)
        f = saturation(A[i])
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

# Layer occupancy
occupancy(l; f=isone) = length(findall(f, l)) / length(l.grid)

# NDI
ndi(x, y) = (x - y) / (x + y)

# Efficiency integral & difference
efficiency_integral(n, k=10_000) = ((n * log(n) - n * log(n + k) + k))

# Efficiency
function efficiency(x, y; k=10_000, rmse=false, kw...)
    n, _rmse = efficiency_gridsearch(x, y; rmse=true, kw...)
    ei = efficiency_integral(exp(n), k)
    if rmse
        return (; ei=ei, rmse=_rmse)
    else
        return ei
    end
end

# Efficiency difference
function efficiency_difference(n, n2; k=10_000)
    n = exp(n)
    n2 = exp(n2)
    return efficiency_integral(n, k) - efficiency_integral(n2, k)
end

# Compare efficiencies within simulations
function comparewithin(
    effs_combined, set; to=nothing, labels=Dict(set .=> set), f=(x, y) -> -(x, y)
)
    # Make sure all labels are defined
    all_labels = Dict(s in keys(labels) ? s => labels[s] : s => s for s in set)
    # Select variables in set & prepare for comparison
    within_combined = @chain effs_combined begin
        @rsubset(:variable in set)
        @rtransform(:variable = all_labels[:variable])
        unstack(:variable, :eff)
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
    for (c1, c2) in comps
        l1 = all_labels[c1]
        l2 = all_labels[c2]
        complabel = "Δ$(l1)_$(l2)"
        @rtransform!(within_combined, $complabel = f($(l1), $(l2)))
    end
    # Select only the variables from the comparison
    within_combined = @chain within_combined begin
        select(:sim, :set, :occ, r"Δ")
        stack(r"Δ")
        dropmissing()
    end
    return within_combined
end

# Flip some comparison values
function flipthatcomp!(df, toflip)
    for comp in toflip
        # Arrange new comparison
        v1, v2 = split(replace(comp, "Δ" => ""), "_")
        newcomp = "Δ$(v2)_$(v1)"
        # Get current indices
        inds = findall(==(comp), df.variable)
        # Update
        new = @view df[inds, :]
        @rtransform!(new, :variable = newcomp)
        @rtransform!(new, :value = -(:value))
    end
    return df
end
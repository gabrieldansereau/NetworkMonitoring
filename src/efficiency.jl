# Saturation equation
saturation(a) = (x) -> x ./ (exp(a) .+ x)

# Efficiency grid search
function efficiency(x, y; A=LinRange(-5.0, 15.0, 10_000))
    err = zeros(length(A))
    for i in eachindex(A)
        f = saturation(A[i])
        err[i] = sqrt(sum((y .- f(x)) .^ 2.0))
    end
    return A[last(findmin(err))]
end

# Layer occupancy
occupancy(l; f=isone) = length(findall(f, l)) / length(l)

# NDI
ndi(x, y) = (x - y) / (x + y)

# Efficiency difference
function efficiency_difference(n, n2; k=10_000)
    n = exp(n)
    n2 = exp(n2)
    return ((n * log(n) - n * log(n + k) + k) - (n2 * log(n2) - n2 * log(n2 + k) + k))
end

# Compare efficiencies within simulations
function comparewithin(effs_combined, set; labels=Dict(set .=> set), f=(x, y) -> -(x, y))
    # Select variables in set & prepare for comparison
    within_combined = @chain effs_combined begin
        @rsubset(:variable in set)
        @rtransform(:variable = :variable in keys(labels) ? labels[:variable] : :variable)
        unstack(:variable, :eff)
    end
    # Compare unique combinations
    labs = [s in keys(labels) ? labels[s] : s for s in set]
    for (l1, l2) in collect(combinations(labs, 2))
        complabel = "Δ$(l1)_$(l2)"
        @rtransform!(within_combined, $complabel = f($(l1), $(l2)))
    end
    # Select comparison variables only
    within_combined = @chain within_combined begin
        select(:sim, :set, :occ, r"Δ")
        stack(r"Δ")
        dropmissing()
    end
    return within_combined
end
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

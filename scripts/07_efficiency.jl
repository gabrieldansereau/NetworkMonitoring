using DrWatson
@quickactivate :NetworkMonitoring

update_theme!(; CairoMakie=(; px_per_unit=2.0))

## Load simulations

# Load summarized simulation results
sims_samplers = CSV.read(datadir("sims_samplers.csv"), DataFrame)
sims_optimized = CSV.read(datadir("sims_optimized.csv"), DataFrame)

# Select UncertaintySampling only as example
sims_set = @rsubset(sims_samplers, :sampler == "UncertaintySampling")

# Visualize
begin
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites", ylabel="Proportion")
    for r in unique(sims_set.sim)
        _sim = @rsubset(sims_set, :sim == r)
        scatter!(ax, _sim.nbon, _sim.med; color=Makie.wong_colors()[2])
    end
    fig
end

## Efficiency

# Simplify data
X = Matrix(select(sims_set, [:sim, :nbon, :med]))

# Saturation equation
saturation(a) = (x) -> x ./ (a .+ x)

# Efficiency grid search
function efficiency(x, y; A=LinRange(-12.0, 12.0, 500))
    err = zeros(length(A))
    for i in eachindex(A)
        f = saturation(exp(A[i]))
        err[i] = sqrt(sum((y .- f(x)) .^ 2.0))
    end
    return exp(A[last(findmin(err))])
end

# Visualize
begin
    f = Figure()
    ax = Axis(f[1, 1])
    for iter in unique(X[:, 1])
        X₁ = X[X[:, 1] .== iter, :][:, 2:3]
        x = X₁[:, 1]
        y = X₁[:, 2]
        eff = efficiency(x, y)
        scatter!(ax, x, y; color=:lightgrey)
        lines!(ax, x, saturation(eff)(x); color=:black, linestyle=:dash)
    end
    f
end

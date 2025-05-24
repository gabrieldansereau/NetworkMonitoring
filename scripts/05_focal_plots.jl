using DrWatson
@quickactivate :NetworkMonitoring

update_theme!(; CairoMakie=(; px_per_unit=2.0))

## Load focal species results

# Load all results
monitored_types = CSV.read(datadir("monitored_types.csv"), DataFrame)
monitored_types2 = CSV.read(datadir("monitored_types2.csv"), DataFrame)
monitored_spp = CSV.read(datadir("monitored_spp.csv"), DataFrame)
monitored_samplers = CSV.read(datadir("monitored_samplers.csv"), DataFrame)
monitored_optimized = CSV.read(datadir("monitored_optimized.csv"), DataFrame)

# Get species with highest degree
sp = monitored_types.sp[1]
deg = maximum(monitored_types.deg)

## Define labels and colors for all plots

# Colors
cols = Dict{Any,Any}(
    # Interaction types
    "possible" => Makie.wong_colors()[2],
    "realized" => Makie.wong_colors()[3],
    "detected" => Makie.wong_colors()[4],
    # Samplers
    "UncertaintySampling" => Makie.wong_colors()[2],
    "WeightedBalancedAcceptance" => Makie.wong_colors()[3],
    "SimpleRandom" => Makie.wong_colors()[1],
    # Layers
    "Focal species range" => Makie.wong_colors()[2],
    "Species richness" => Makie.wong_colors()[4],
    "Realized interactions" => Makie.wong_colors()[5],
)
for (sp, col) in zip(unique(monitored_spp.sp), [Makie.wong_colors()[[2, 6, 7]]..., :black])
    cols[sp] = col
end
cols

## Monitored types

begin
    res = monitored_types
    vars = unique(res.var)
    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions", xticks=0:25:100
    )
    for v in vars
        b = filter(:var => ==(v), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(b.nbon, b.med; label=v, color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rc, merge=true)
    fig
end
save(plotsdir("focal_types.png"), fig)

# Visualize result
begin
    res = monitored_types2
    vars = unique(res.var)
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel="Sites in BON",
        ylabel="Monitored interactions",
        xticks=0:2000:10_000,
    )
    for v in vars
        b = filter(:var => ==(v), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=v, color=cols[v])
        lines!(b.nbon, b.med; label=v, color=cols[v])
        hlines!(unique(b.deg); linestyle=:dash, color=cols[v])
    end
    hlines!(ax, maximum(res.deg); linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:rb, merge=true)
    fig
end
save(plotsdir("focal_types2.png"), fig)

## Repeat with 4 species with different degrees

# Visualize result
begin
    res = monitored_spp
    spp = unique(res.sp)
    fig = Figure(; size=(600, 600))
    ax1 = Axis(fig[1, 1]; ylabel="Monitored interactions", xticks=0:20:100)
    ax2 = Axis(
        fig[2, 1];
        ylabel="Proportion monitored",
        xlabel="Number of sites in BON",
        xticks=0:20:100,
        yticks=0:0.2:1.0,
    )
    for (i, sp) in enumerate(spp)
        b = filter(:sp => ==(sp), res)
        # Monitored int
        band!(ax1, b.nbon, b.low, b.upp; alpha=0.4, label=sp)
        lines!(ax1, b.nbon, b.med; label=sp)
        hlines!(ax1, unique(b.deg); linestyle=:dash, alpha=0.5)
        # Proportion
        band!(ax2, b.nbon, b.low ./ b.deg, b.upp ./ b.deg; alpha=0.4, label=sp)
        lines!(ax2, b.nbon, b.med ./ b.deg; label=sp)
    end
    fig[:, end + 1] = Legend(fig, ax1, "Species"; framevisible=false, merge=true)
    fig
end
save(plotsdir("focal_spp.png"), fig)

## Explore variations with different sampler

# Plot
begin
    res = monitored_samplers
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=s, color=cols[s])
        lines!(b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, deg; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend(; position=:lt, merge=true)
    fig
end
save(plotsdir("focal_samplers.png"), fig)

## Optimized sampling

# Plot
begin
    res = monitored_optimized
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sites in BON", ylabel="Monitored interactions")
    for s in unique(res.sampler)
        b = filter(:sampler => ==(s), res)
        band!(b.nbon, b.low, b.upp; alpha=0.4, label=s, color=cols[s])
        lines!(b.nbon, b.med; label=s, color=cols[s])
    end
    hlines!(ax, deg; linestyle=:dash, alpha=0.5, color=:grey, label="metaweb")
    axislegend("Optimized on"; position=:lt, merge=true)
    fig
end
save(plotsdir("focal_optimized.png"), fig)

### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ f8c83bf9-7560-4890-8f10-2dd03cd00b7a
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    using Revise
    using DrWatson
    using Printf
    using NetworkMonitoring
    update_theme!(; CairoMakie=(; px_per_unit=2.0))
end

# ╔═╡ bec626f6-c174-4988-97be-bbb992d16890
md"""
2026-03-07

Investigating number of replicates and change in efficiency curve
"""

# ╔═╡ e0098333-f8b5-4413-ad63-343c8c4d2fc8
md"""
## Load focal species results
"""

# ╔═╡ 9da4af9b-07f2-4bf0-8df6-2fefeb8778f8
begin
    ## Load focal species results

    # Use job id to vary parameters
    id = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    idp = lpad(id, 3, "0")

    # Load & summarize test results
    monitored_test_all = CSV.read(datadir("monitored_test.csv"), DataFrame)
    monitored_test = summarize_focal(monitored_test_all; id=id)
end

# ╔═╡ 1ddf24c4-d555-4cff-8381-84359ebc804f
md"""
## Range estimation
"""

# ╔═╡ f3ccd2fb-88c6-40ae-99cb-9b3872f54138
md"""
### Define functions
"""

# ╔═╡ 61719ed0-9ed8-45a1-b76a-60e673f1b0a2
# Load layers used for optimization
function load_layers(idp)
    errors = -0.2:0.05:0.2
    estimated_ranges = Dict()
    for (i, e) in enumerate(errors)
        estimated_ranges[e] = SDT.SDMLayer(
            datadir("focal_array", "layer_range_estimations-$idp.tiff"); bandnumber=i
        )
    end
    return estimated_ranges
end

# ╔═╡ 31614fca-b5b1-4c93-b696-b0ad7a783cc2
# The set might differ, so just in case
function fix_set(estimated_ranges)
    set = ["Over-0.2", "True-0.0", "Under-0.2"]
    for (s, v) in zip(set, [-0.2, 0.0, 0.2])
        estimated_ranges[s] = estimated_ranges[v]
    end
    return estimated_ranges
end

# ╔═╡ f7428070-4e29-4522-8eda-76f4627303d2
md"""
### Load simulation
"""

# ╔═╡ dbf0dca1-9eb6-47cf-a750-d3635b87c0ab
estimated_ranges = load_layers(idp)

# ╔═╡ 2f0a2245-ca81-4b81-9cd4-eb6ce86bc847
fix_set(estimated_ranges)

# ╔═╡ 99bc889c-068b-43f7-869a-1666a5d07d57
md"""
### Check results
"""

# ╔═╡ 25055648-5d6b-4e14-be72-5a0b40364f60
OUTDIR = "dev";

# ╔═╡ 87db2694-a821-424e-98e9-342624dbffb6
# Load & summarize results
function load_sim(exp; idp=idp)
    monitored_estimations_all = CSV.read(
        datadir(OUTDIR, "monitored_estimations-$idp-$exp.csv"), DataFrame
    )
    id = parse(Int, idp)
    return monitored_estimations = summarize_focal(monitored_estimations_all; id=id)
end

# ╔═╡ 4bfd447b-b234-4ee7-83a5-abd57d86d224
NREP = 50;

# ╔═╡ 08db2a32-45fd-4839-bb4d-e06483d63713
STEP = 05;

# ╔═╡ cc54a065-4a01-422f-bc34-c93816a37bd8
EXP = "$(lpad(NREP,3,"0"))-$(lpad(STEP,2,"0"))";

# ╔═╡ 3f64b170-4b62-4786-ac62-ddfe6714d96f
begin
    monitored_estimations = load_sim(EXP)
    first(monitored_estimations, 3)
end

# ╔═╡ 8c9934a6-77ea-43a9-8d20-35e3b40d1503
set = ["Over-0.2", "True-0.0", "Under-0.2"];

# ╔═╡ a3bf9f98-f377-4562-8b8e-dbd81819c6ba
# Generate BON examples
begin
    Random.seed!(33)
    bons = Dict()
    for e in set
        bons[e] = BON.sample(BON.BalancedAcceptance(100), estimated_ranges[e])
    end
    bons
end

# ╔═╡ e9c58db3-8a2b-4bd6-b578-95e6412f59f4
md"""
## Check multiple options
"""

# ╔═╡ 172082a8-acb4-40d1-92b0-65f849475cd8
md"""
#### STEP = 5
"""

# ╔═╡ 7401ebb7-65ee-46f5-9f66-c5b0b98a5c9b
md"""
#### STEP = 10
"""

# ╔═╡ 6926fbac-dfbd-4b88-916b-4f8eb4be2565
md"""
## Investigate error and efficiency variation with number of replicates and step
"""

# ╔═╡ 74d920c3-c07e-4513-b289-f40217030efb
rmse_df = DataFrame();

# ╔═╡ 48ef09f2-279e-49ef-a794-720ddf65303c
# Plot
function plot_sim(
    monitored_estimations, estimated_ranges, EXP; show_eff=true, scatter=false, push=true
)
    return fig_estimation = let
        set = ["Over-0.2", "True-0.0", "Under-0.2"]
        var = :layer
        res = filter(var => in(set), monitored_estimations)
        vals = unique(res[:, var])

        range_over = estimated_ranges[set[1]]
        range_true = estimated_ranges[set[2]]
        range_under = estimated_ranges[set[3]]

        if !(@isdefined colours)
            colours = Dict()
        end
        for (i, s) in enumerate(set)
            colours[s] = Makie.wong_colors()[i]
        end

        # Create figure
        fig = Figure()
        # Create layouts
        ga = GridLayout(fig[:, 1:3])
        gb = GridLayout(fig[:, end + 1])
        # Create axes
        ax = Axis(
            ga[1, 1];
            xlabel="Sites in BON",
            ylabel="Monitored interactions",
            xticks=0:100:500,
        )
        ax1 = Axis(
            gb[1, 1];
            aspect=1,
            yaxisposition=:right,
            ylabelrotation=1.5pi,
            ylabelsize=10,
        )
        ax2 = Axis(
            gb[2, 1];
            aspect=1,
            yaxisposition=:right,
            ylabelrotation=1.5pi,
            ylabelsize=10,
        )
        ax3 = Axis(
            gb[3, 1];
            aspect=1,
            yaxisposition=:right,
            ylabelrotation=1.5pi,
            ylabelsize=10,
        )
        # Remove decorations for heatmaps
        hidedecorations!(ax1; label=false)
        hidedecorations!(ax2; label=false)
        hidedecorations!(ax3; label=false)

        # Sampling results
        nrep, step = split(EXP, "-")
        for v in vals
            b = filter(var => ==(v), res)
            eff_a = efficiency_gridsearch(b.nbon, b.med; f=exp)
            eff, rmse = efficiency(b.nbon, b.med; rmse=true, f=exp)
            eff_low, _ = efficiency(b.nbon, b.med .- rmse; rmse=true, f=exp)
            eff_upp, _ = efficiency(b.nbon, b.med .+ rmse; rmse=true, f=exp)
            lab = "$v, eff=$(@sprintf("%.3f", eff)), rmse=$(@sprintf("%.5f", rmse))"
            if push
                push!(
                    rmse_df,
                    (
                        id=parse(Int, idp),
                        nrep=nrep,
                        step=step,
                        value=v,
                        eff=eff,
                        log_a=log(eff_a),
                        rmse=rmse,
                        low=eff_low,
                        upp=eff_upp,
                        range=eff_upp - eff_low,
                    ),
                )
            end
            band!(ax, b.nbon, b.low, b.upp; alpha=0.4, color=colours[v], label=lab)
            if show_eff
                lines!(
                    ax,
                    b.nbon,
                    saturation(eff_a)(b.nbon);
                    color=:black,
                    linestyle=:dash,
                    linewidth=1.5,
                    alpha=0.7,
                )
            end
            if scatter
                scatter!(ax, b.nbon, b.med; color=colours[v], label=lab, markersize=5.0)
            end
            lines!(ax, b.nbon, b.med; label=lab)
        end
        hlines!(ax, [1.0]; linestyle=:dash, alpha=0.5, color=:grey)
        Legend(ga[2, 1], ax; orientation=:horizontal, merge=true, nbanks=3)
        # Heatmaps & BON example
        heatmap!(ax1, range_over; colormap=:greys, alpha=0.5)
        heatmap!(ax1, range_true; colormap=:viridis)
        heatmap!(ax2, range_true; colormap=:viridis)
        heatmap!(ax3, range_true; colormap=:viridis, alpha=0.5)
        heatmap!(ax3, range_under; colormap=:viridis)
        for (a, v) in zip([ax1, ax2, ax3], vals)
            scatter!(a, coordinates(bons[v]); markersize=5, strokewidth=0.5, color=colours[v])
            a.ylabel = v
        end

        # Subpanel labels
        Label(
            ga[1, :, Top()],
            "Simulation $idp with NREP = $nrep, STEP=$step";
            padding=(0, 0, 5, 0),
            font=:bold,
        )
        Label(gb[1, :, Top()], "BON examples"; padding=(0, 0, 5, 0), font=:bold)
        # Show figure
        save(plotsdir("_xtras", "focal_ranges-$EXP.png"), fig)
        fig
    end
end

# ╔═╡ e719aef9-9309-497a-a07a-a6747e07265c
function load_and_plot(EXP, idp)
    monitored_estimations = load_sim(EXP)
    return plot_sim(monitored_estimations, estimated_ranges, EXP)
end

# ╔═╡ 652d6511-8a73-4e47-8388-5001d6195dac
load_and_plot("100-05", "01")

# ╔═╡ ceebfdff-d390-41cd-9588-2a23a0dfe6a9
load_and_plot("050-05", "01")

# ╔═╡ ee411f14-c9d8-4df6-9bcc-d04a345ceedc
load_and_plot("030-05", "01")

# ╔═╡ e3166c34-ca92-45c0-8b92-151aea6918b9
load_and_plot("020-05", "01")

# ╔═╡ 0dfdfc97-c9f4-4da6-8d4c-8cb4170851d1
load_and_plot("050-10", "01")

# ╔═╡ 672ceea3-6785-4b63-9f44-146246d2cd92
load_and_plot("030-10", "01")

# ╔═╡ bea2686e-0a1a-4413-a739-e04fb5943d75
load_and_plot("020-10", "01")

# ╔═╡ 9345d8c8-0f1d-4fa3-8c3b-d5c182a07676
plot_sim(monitored_estimations, estimated_ranges, EXP; push=false)

# ╔═╡ 6f56f00c-3bf6-4a45-a784-c922dd672124
sort(@rsubset(rmse_df, :value == "True-0.0"), :rmse)

# ╔═╡ 503249ac-1380-44e9-a21b-ce82b2253c25
sort(@rsubset(rmse_df, :value == "Under-0.2"), :rmse)

# ╔═╡ 60f91fc3-9809-40d4-887f-0e0723d19bb2
sort(@rsubset(rmse_df, :value == "Over-0.2"), :rmse)

# ╔═╡ 6b00bc1d-02cf-41da-83f6-a1fd3a0904d5
CSV.write(datadir("rmse_investigation.csv"), rmse_df);

# ╔═╡ Cell order:
# ╟─bec626f6-c174-4988-97be-bbb992d16890
# ╠═f8c83bf9-7560-4890-8f10-2dd03cd00b7a
# ╟─e0098333-f8b5-4413-ad63-343c8c4d2fc8
# ╟─9da4af9b-07f2-4bf0-8df6-2fefeb8778f8
# ╟─1ddf24c4-d555-4cff-8381-84359ebc804f
# ╟─f3ccd2fb-88c6-40ae-99cb-9b3872f54138
# ╟─87db2694-a821-424e-98e9-342624dbffb6
# ╟─61719ed0-9ed8-45a1-b76a-60e673f1b0a2
# ╟─31614fca-b5b1-4c93-b696-b0ad7a783cc2
# ╠═48ef09f2-279e-49ef-a794-720ddf65303c
# ╟─e719aef9-9309-497a-a07a-a6747e07265c
# ╟─f7428070-4e29-4522-8eda-76f4627303d2
# ╠═3f64b170-4b62-4786-ac62-ddfe6714d96f
# ╠═dbf0dca1-9eb6-47cf-a750-d3635b87c0ab
# ╠═2f0a2245-ca81-4b81-9cd4-eb6ce86bc847
# ╟─a3bf9f98-f377-4562-8b8e-dbd81819c6ba
# ╟─99bc889c-068b-43f7-869a-1666a5d07d57
# ╠═25055648-5d6b-4e14-be72-5a0b40364f60
# ╠═4bfd447b-b234-4ee7-83a5-abd57d86d224
# ╠═08db2a32-45fd-4839-bb4d-e06483d63713
# ╠═cc54a065-4a01-422f-bc34-c93816a37bd8
# ╠═8c9934a6-77ea-43a9-8d20-35e3b40d1503
# ╠═9345d8c8-0f1d-4fa3-8c3b-d5c182a07676
# ╟─e9c58db3-8a2b-4bd6-b578-95e6412f59f4
# ╟─172082a8-acb4-40d1-92b0-65f849475cd8
# ╟─652d6511-8a73-4e47-8388-5001d6195dac
# ╟─ceebfdff-d390-41cd-9588-2a23a0dfe6a9
# ╟─ee411f14-c9d8-4df6-9bcc-d04a345ceedc
# ╟─e3166c34-ca92-45c0-8b92-151aea6918b9
# ╟─7401ebb7-65ee-46f5-9f66-c5b0b98a5c9b
# ╟─0dfdfc97-c9f4-4da6-8d4c-8cb4170851d1
# ╟─672ceea3-6785-4b63-9f44-146246d2cd92
# ╟─bea2686e-0a1a-4413-a739-e04fb5943d75
# ╟─6926fbac-dfbd-4b88-916b-4f8eb4be2565
# ╠═74d920c3-c07e-4513-b289-f40217030efb
# ╠═6f56f00c-3bf6-4a45-a784-c922dd672124
# ╠═503249ac-1380-44e9-a21b-ce82b2253c25
# ╠═60f91fc3-9809-40d4-887f-0e0723d19bb2
# ╠═6b00bc1d-02cf-41da-83f6-a1fd3a0904d5

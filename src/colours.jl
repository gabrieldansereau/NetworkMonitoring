# Define colors for all plots
colours = Dict{Any,Any}(
    # Interaction types
    "possible" => Makie.wong_colors()[2],
    "realized" => Makie.wong_colors()[3],
    "detected" => Makie.wong_colors()[4],
    # Samplers
    "Uncertainty Sampling" => Makie.wong_colors()[2],
    "Weighted Balanced Acceptance" => Makie.wong_colors()[3],
    "Simple Random" => Makie.wong_colors()[1],
    "Balanced Acceptance" => Makie.wong_colors()[1],
    "Balanced Acceptance Mask" => :grey,
    "Simple Random Mask" => :turquoise,
    # Layers
    "Focal species range" => Makie.wong_colors()[2],
    "Species richness" => Makie.wong_colors()[4],
    "Realized interactions" => Makie.wong_colors()[5],
    "Probabilistic range" => Makie.wong_colors()[6],
)

# AlgebraOfGraphics needs a different format
colourpal = [k => v for (k, v) in colours]
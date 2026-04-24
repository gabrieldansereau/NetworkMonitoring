function _median_confint(x; α=0.05)
    z = quantile(Distributions.Normal(0.0, 1.0), 1 - α / 2)
    n = length(x)
    L = ceil(Int, 0.5 * n - z * sqrt(0.25 * n))
    U = ceil(Int, 0.5 * n + z * sqrt(0.25 * n))
    return (low=sort(x)[L], upp=sort(x)[U])
end

function summarize_focal(df; id=0, confint=false, α=0.05)
    notcols = [:rep, :monitored]
    cols = setdiff(Symbol.(names(df)), notcols)
    monitored = @chain df begin
        groupby(Not(notcols))
        @combine(
            :low = quantile(:monitored, 0.05),
            :med = median(:monitored),
            :upp = quantile(:monitored, 0.95),
            :deg = maximum(:deg)
        )
        @rtransform(:low = :low / :deg, :med = :med / :deg, :upp = :upp / :deg,)
        @select(:sim = id, All())
    end
    if confint
        monitored_confint = @chain df begin
            groupby(Not(notcols))
            @combine(
                :deg = maximum(:deg),
                :confint_low = Ref(:monitored),
                :confint_upp = Ref(:monitored),
            )
            @rtransform(
                :confint_low = _median_confint(:confint_low; α=α).low,
                :confint_upp = _median_confint(:confint_upp; α=α).upp,
            )
            @rtransform(
                :confint_low = :confint_low / :deg, :confint_upp = :confint_upp / :deg,
            )
            @select(Not(:deg))
        end
        leftjoin!(
            monitored,
            monitored_confint;
            on=intersect(cols, Symbol.(names(monitored_confint))),
        )
        ordered = [:sim, :layer, :offset, :nbon, :deg, :degmax, :low, :med, :upp]
        select!(monitored, filter(in(names(monitored)), ordered), r"^confint", All())
    end
    monitored.sampler =
        replace.(
            monitored.sampler,
            "UncertaintySampling" => "Uncertainty Sampling",
            "WeightedBalancedAcceptance" => "Weighted Balanced Acceptance",
            "BalancedAcceptanceMask" => "Balanced Mask",
            "BalancedAcceptance" => "Balanced Acceptance",
            "SimpleRandomMask" => "Simple Random Mask",
            "SimpleRandom" => "Simple Random",
        )
    return monitored
end
